from cereal import car
from common.numpy_fast import clip, interp
from common.params import Params
from common.realtime import DT_CTRL
from selfdrive.controls.lib.drive_helpers import CONTROL_N, apply_deadzone
from selfdrive.controls.lib.pid import PIDController
from selfdrive.modeld.constants import T_IDXS

LongCtrlState = car.CarControl.Actuators.LongControlState
ButtonType = car.CarState.ButtonEvent.Type

# 정차 래치. 완전히 선 뒤에는 플래너가 한 프레임 튀어도 stopping 을 유지한다.
STANDSTILL_LATCH_SPEED = 0.05
# 출발은 플래너가 이만큼 "연속으로" 요구해야 인정한다(DT_CTRL 0.01 s 기준 0.5 s).
# 이전 값 10 프레임(0.1 s)은 리드가 잠깐 끊기거나 정지거리 오차를 메우려는
# 플래너 요구만으로도 래치가 풀려, 앞차가 서 있는데 전진하는 원인이 되었다.
STANDSTILL_RELEASE_FRAMES = 50

# 정차 중 앞차 출발 판정.
# 앞차가 보이는 동안에는 플래너 요구만으로 출발하지 않고, 앞차가 실제로
# 멀어지고 있을 때만 출발 요구를 인정한다.
LEAD_DEPART_SPEED_RADAR = 0.5    # m/s, 레이더 리드 속도
LEAD_DEPART_SPEED_VISION = 1.0   # m/s, 비전 리드는 속도 잡음이 커서 더 높게
LEAD_DEPART_DIST_RADAR = 1.0     # m, 정차 기준거리 대비 증가량
LEAD_DEPART_DIST_VISION = 2.0    # m, 근거리 비전 거리는 몇 m씩 튀므로 더 크게
# 정차할 때 리드가 있었는데 사라진 경우, 이만큼 "계속" 안 보여야 리드가 떠난
# 것으로 본다(2.0 s). 근접 정차에서 레이더·모델이 리드를 잠깐 놓치는 것만으로
# 출발하지 않게 한다. 운전자 가속페달·RES 는 이와 무관하게 즉시 출발한다.
LEAD_LOST_RELEASE_FRAMES = 200

def long_control_state_trans(CP, active, long_control_state, v_ego, v_target,
                             v_target_1sec, brake_pressed, cruise_standstill,
                             soft_hold, a_target_now, starting_state,
                             standstill_latched=False, standstill_release=False):
  # apilot-c2 stopping transition: keep PID braking while the planned
  # acceleration is still strong, then hand over to the stopping ramp.
  # Match aPilot C2: CarState suppresses the stock SCC standstill flag during
  # openpilot longitudinal control and retains it for stock ACC.
  cruise_standstill = cruise_standstill and not CP.enableGasInterceptor
  accelerating = v_target_1sec > (v_target + 0.01)
  planned_stop = (v_target < CP.vEgoStopping and
                  v_target_1sec < CP.vEgoStopping and
                  not accelerating)
  stay_stopped = (v_ego < CP.vEgoStopping and
                  (brake_pressed or cruise_standstill))
  stopping_condition = planned_stop or stay_stopped

  starting_condition = (v_target_1sec > CP.vEgoStarting and
                        accelerating and
                        not cruise_standstill and
                        not brake_pressed)
  started_condition = v_ego > CP.vEgoStarting

  if not active:
    long_control_state = LongCtrlState.off

  else:
    if long_control_state in (LongCtrlState.off, LongCtrlState.pid):
      long_control_state = LongCtrlState.pid
      if stopping_condition and a_target_now > -1.0:
        long_control_state = LongCtrlState.stopping

    elif long_control_state == LongCtrlState.stopping:
      if starting_condition or standstill_release:
        long_control_state = LongCtrlState.starting if starting_state else LongCtrlState.pid

    elif long_control_state == LongCtrlState.starting:
      if stopping_condition:
        long_control_state = LongCtrlState.stopping
      elif started_condition:
        long_control_state = LongCtrlState.pid

    # 완전히 선 뒤에는 플래너 한 프레임의 잡음으로 SCC StopReq 를 놓지 않는다.
    # 확정된 출발(standstill_release)만 stopping 을 벗어날 수 있다.
    if standstill_latched and not standstill_release:
      long_control_state = LongCtrlState.stopping

    if soft_hold:
      long_control_state = LongCtrlState.stopping

  return long_control_state, planned_stop


class LongControl:
  def __init__(self, CP):
    self.CP = CP
    self.long_control_state = LongCtrlState.off
    self.pid = PIDController((CP.longitudinalTuning.kpBP, CP.longitudinalTuning.kpV),
                             (CP.longitudinalTuning.kiBP, CP.longitudinalTuning.kiV),
                             k_f=CP.longitudinalTuning.kf, rate=1 / DT_CTRL)
    self.params = Params()
    self.read_param_count = 0
    self.stop_accel = CP.stopAccel
    self.long_coast_band = 0.0
    self.v_pid = 0.0
    self.last_output_accel = 0.0
    self.starting_accel = 0.0
    self.starting_ramp_rate = 2.0
    self.standstill_hold_accel = -1.1
    self.standstill_hold_rate = 1.2
    self.standstill_latched = False
    self.standstill_frames = 0
    # 정차 중 앞차 추적: 래치 중 리드를 본 적 있는지, 출발 판정 기준거리,
    # 리드가 연속으로 안 보인 프레임 수.
    self.standstill_lead_seen = False
    self.standstill_lead_ref_d = None
    self.standstill_lead_lost_frames = 0
    # 정차 중 실제로 걸고 있던 제동값. stopping 에 다시 들어올 때 0 부터 램프를
    # 다시 타지 않게 해서, 상태가 잠깐 흔들려도 제동이 약해지지 않는다.
    self.standstill_hold_memory = None
    self._update_pid_gains()
    # Read launch control immediately so StartAccelApply=0 disables the
    # starting state from the first control cycle.
    self._update_start_accel()
    self._update_stop_accel()
    self._update_stopping_decel_rate()
    self._update_standstill_hold()

    # apilot-c2 uses two actuator-delay predictions and selects the more
    # conservative target. Derive safe defaults around the configured delay.
    delay = float(clip(CP.longitudinalActuatorDelay, 0.1, 1.0))
    schema_lower = CP.longitudinalActuatorDelayLowerBound
    schema_upper = CP.longitudinalActuatorDelayUpperBound
    self.actuator_delay_lower = float(clip(schema_lower if schema_lower > 0.0 else delay - 0.1, 0.1, 0.99))
    self.actuator_delay_upper = float(clip(schema_upper if schema_upper > 0.0 else delay + 0.1,
                                           self.actuator_delay_lower, 1.0))
    self._update_actuator_delays()

  def _update_start_accel(self):
    start_raw = self.params.get_int("StartAccelApply", 0)
    self.start_accel_apply = float(clip(start_raw * 0.01, 0.0, 1.0))
    self.start_accel = float(clip(2.0 * self.start_accel_apply, 0.0, 2.0))
    self.starting_state = start_raw > 0

  def _update_stop_accel(self):
    stop_raw = self.params.get("StopAccelApply", encoding="utf8")
    if stop_raw is not None:
      try:
        stop_accel_apply = float(clip(int(stop_raw) * 0.01, 0.0, 1.0))
      except (TypeError, ValueError):
        stop_accel_apply = 0.3
      self.stop_accel = -2.0 * stop_accel_apply
    else:
      # Preserve an existing StoppingAccel value until StopAccelApply is
      # changed in the UI. With neither value set, use the car default.
      legacy_stop_accel = self.params.get_float("StoppingAccel") * 0.01
      self.stop_accel = legacy_stop_accel if legacy_stop_accel < 0.0 else self.CP.stopAccel

  def _update_stopping_decel_rate(self):
    rate_raw = self.params.get_int("StoppingDecelRate")
    rate = rate_raw * 0.01 if rate_raw > 0 else self.CP.stoppingDecelRate
    self.stopping_decel_rate = float(clip(rate, 0.2, 2.0))

  def _reset_standstill_lead(self):
    self.standstill_lead_seen = False
    self.standstill_lead_ref_d = None
    self.standstill_lead_lost_frames = 0

  def _standstill_lead_allows_start(self, radar_state):
    """정차 중 앞차 상태로 보아 플래너의 출발 요구를 받아들여도 되는지.

    - 정차한 뒤 리드를 한 번도 못 봤으면(신호 대기 등) 플래너를 따른다.
    - 리드가 보이면 실제로 멀어질 때만 허용한다. 정지거리를 메우려는
      "간격 채우기" 요구는 여기서 막힌다.
    - 리드를 봤는데 사라졌으면 LEAD_LOST_RELEASE_FRAMES 동안 계속 안 보일
      때만 허용한다. 근접 정차에서 리드가 잠깐 끊기는 것을 걸러 낸다.
    """
    lead = radar_state.leadOne if radar_state is not None else None
    if lead is not None and lead.status:
      if not self.standstill_lead_seen or self.standstill_lead_ref_d is None:
        # 처음 잡았거나 끊겼다가 다시 잡은 리드는 지금 거리를 기준으로 삼는다.
        self.standstill_lead_ref_d = float(lead.dRel)
      self.standstill_lead_seen = True
      self.standstill_lead_lost_frames = 0
      if lead.radar:
        depart_speed, depart_dist = LEAD_DEPART_SPEED_RADAR, LEAD_DEPART_DIST_RADAR
      else:
        depart_speed, depart_dist = LEAD_DEPART_SPEED_VISION, LEAD_DEPART_DIST_VISION
      return (lead.vLead > depart_speed or
              lead.dRel - self.standstill_lead_ref_d > depart_dist)

    self.standstill_lead_ref_d = None
    if not self.standstill_lead_seen:
      return True
    self.standstill_lead_lost_frames += 1
    return self.standstill_lead_lost_frames >= LEAD_LOST_RELEASE_FRAMES

  def _update_standstill_latch(self, active, CS, v_target, v_target_1sec, soft_hold, radar_state=None):
    """완전 정차를 래치하고, 확정된 출발 요구에서만 푼다.

    플래너 요구만 보던 예전 방식은 근접 정차에서 리드가 잠깐 끊기거나,
    목표 정지거리보다 일찍 섰을 때 간격을 메우려는 요구만으로 래치가 풀려
    앞차가 서 있는데 전진했다. 이제 앞차가 실제로 출발했는지(또는 충분히
    오래 사라졌는지)를 함께 확인하고, 그 상태가 0.5 초 이어져야 출발한다.
    """
    if not active:
      self.standstill_latched = False
      self.standstill_frames = 0
      self.standstill_hold_memory = None
      self._reset_standstill_lead()
      return False, False

    if self.long_control_state == LongCtrlState.stopping and CS.vEgo < STANDSTILL_LATCH_SPEED:
      self.standstill_latched = True

    if not self.standstill_latched:
      self.standstill_frames = 0
      self._reset_standstill_lead()
      return False, False

    planner_start = (v_target_1sec > self.CP.vEgoStarting and
                     v_target_1sec > v_target + 0.01 and
                     not CS.brakePressed and
                     not CS.cruiseState.standstill)
    # 리드 추적은 플래너 요구와 상관없이 매 프레임 갱신해야 사라진 시간을
    # 제대로 센다.
    lead_allows_start = self._standstill_lead_allows_start(radar_state)
    starting_request = planner_start and lead_allows_start
    self.standstill_frames = self.standstill_frames + 1 if starting_request else 0

    resume_pressed = any(
      event.pressed and event.type in (ButtonType.accelCruise, ButtonType.resumeCruise)
      for event in CS.buttonEvents)
    # 운전자 의사(가속페달·RES)는 소프트홀드까지 포함해 언제나 즉시 푼다.
    release = CS.gasPressed or resume_pressed
    if not soft_hold and self.standstill_frames >= STANDSTILL_RELEASE_FRAMES:
      release = True

    if release:
      self.standstill_latched = False
      self.standstill_frames = 0
      self.standstill_hold_memory = None
      self._reset_standstill_lead()
      return False, True

    return True, False

  def _update_standstill_hold(self):
    hold_raw = self.params.get("StandstillHoldApply", encoding="utf8")
    rate_raw = self.params.get("StandstillHoldRate", encoding="utf8")
    try:
      hold_apply = int(hold_raw) if hold_raw is not None else 55
    except (TypeError, ValueError):
      hold_apply = 55
    try:
      hold_rate = int(rate_raw) * 0.01 if rate_raw is not None else 1.2
    except (TypeError, ValueError):
      hold_rate = 1.2

    self.standstill_hold_accel = -2.0 * float(clip(hold_apply * 0.01, 0.1, 1.0))
    self.standstill_hold_rate = float(clip(hold_rate, 0.2, 2.0))

  def _update_actuator_delays(self):
    lower = self.params.get_float("LongitudinalActuatorDelayLowerBound") * 0.01
    upper = self.params.get_float("LongitudinalActuatorDelayUpperBound") * 0.01
    if lower > 0.0:
      self.actuator_delay_lower = float(clip(lower, 0.1, 0.99))
    if upper > 0.0:
      self.actuator_delay_upper = float(clip(upper, self.actuator_delay_lower, 1.0))
    elif self.actuator_delay_upper < self.actuator_delay_lower:
      self.actuator_delay_upper = self.actuator_delay_lower

  def _update_pid_gains(self):
    if len(self.CP.longitudinalTuning.kpBP) != 1 or len(self.CP.longitudinalTuning.kiBP) != 1:
      return

    kp_raw = self.params.get("LongTuningKpV", encoding="utf8")
    ki_raw = self.params.get("LongTuningKiV", encoding="utf8")
    kf_raw = self.params.get("LongTuningKf", encoding="utf8")
    try:
      if kp_raw is not None:
        kp = float(clip(float(kp_raw) * 0.01, 0.0, 2.0))
        self.pid._k_p = (self.CP.longitudinalTuning.kpBP, [kp])
      if ki_raw is not None:
        old_ki = self.pid.k_i
        ki = float(clip(float(ki_raw) * 0.001, 0.0, 2.0))
        self.pid._k_i = (self.CP.longitudinalTuning.kiBP, [ki])
        # The integral is already gain-scaled.  Clear it for Ki=0 and scale it
        # down with a reduced live gain so an old correction cannot keep the
        # accelerator pinned after the setting changes.
        if ki == 0.0:
          self.pid.i = 0.0
        elif old_ki > 0.0 and ki < old_ki:
          self.pid.i *= ki / old_ki
      if kf_raw is not None:
        self.pid.k_f = float(clip(float(kf_raw) * 0.01, 0.0, 2.0))
    except (TypeError, ValueError):
      pass

  def reset(self, v_pid=0.0):
    self.pid.reset()
    self.v_pid = v_pid

  def _read_params(self):
    self.read_param_count += 1
    if self.read_param_count >= 100:
      self.read_param_count = 0
      self._update_stop_accel()
      self._update_stopping_decel_rate()
      self._update_standstill_hold()
      self.long_coast_band = clip(self.params.get_float("LongCoastBand") * 0.01, 0.0, 0.4)
      self._update_actuator_delays()

      self._update_start_accel()

    elif self.read_param_count == 10:
      self._update_pid_gains()

  def update(self, active, CS, long_plan, accel_limits, t_since_plan, soft_hold=False, radar_state=None):
    self._read_params()

    if len(long_plan.speeds) == CONTROL_N:
      speeds = long_plan.speeds
      accels = long_plan.accels
      v_target_now = interp(t_since_plan, T_IDXS[:CONTROL_N], speeds)
      a_target_now = interp(t_since_plan, T_IDXS[:CONTROL_N], accels)
      j_target = long_plan.jerks[0] if len(long_plan.jerks) else 0.0

      # apilot-c2 dual-delay compensation. The lower and upper delay estimates
      # absorb vehicle brake-response variation; use the more braking target.
      v_target_lower = interp(self.actuator_delay_lower + t_since_plan,
                              T_IDXS[:CONTROL_N], speeds)
      v_target_upper = interp(self.actuator_delay_upper + t_since_plan,
                              T_IDXS[:CONTROL_N], speeds)
      a_target_lower = 2.0 * (v_target_lower - v_target_now) / self.actuator_delay_lower - a_target_now
      a_target_upper = 2.0 * (v_target_upper - v_target_now) / self.actuator_delay_upper - a_target_now
      v_target = min(v_target_lower, v_target_upper)
      a_target = min(a_target_lower, a_target_upper)

      v_target_1sec = interp(self.actuator_delay_lower + t_since_plan + 1.0,
                             T_IDXS[:CONTROL_N], speeds)
    else:
      v_target_now = 0.0
      a_target_now = 0.0
      v_target = 0.0
      v_target_1sec = 0.0
      a_target = 0.0
      j_target = 0.0

    self.pid.neg_limit = accel_limits[0]
    self.pid.pos_limit = accel_limits[1]

    prev_long_control_state = self.long_control_state
    standstill_latched, standstill_release = self._update_standstill_latch(
      active, CS, v_target, v_target_1sec, soft_hold, radar_state)
    self.long_control_state, planned_stop = long_control_state_trans(
      self.CP, active, self.long_control_state, CS.vEgo, v_target, v_target_1sec,
      CS.brakePressed, CS.cruiseState.standstill, soft_hold, a_target_now,
      self.starting_state, standstill_latched, standstill_release)

    if self.long_control_state == LongCtrlState.starting and prev_long_control_state != LongCtrlState.starting:
      # Begin every launch from zero drive request instead of jumping straight
      # to StartAccelApply. This removes the pause-then-lurch feeling after the
      # stopping brake is released.
      self.starting_accel = 0.0

    if self.long_control_state == LongCtrlState.off:
      self.reset(CS.vEgo)
      self.starting_accel = 0.0
      output_accel = 0.0

    elif self.long_control_state == LongCtrlState.stopping:
      self.starting_accel = 0.0
      output_accel = self.last_output_accel
      # 정차 중 걸고 있던 제동값을 기억해 두었다면 거기서 이어 간다. 상태가
      # 한 번 흔들렸다고 0 에서 다시 램프를 타면 약 0.9 초 동안 제동이 약해져
      # 차가 조금씩 밀린다("뚝뚝" 풀림).
      if self.standstill_hold_memory is not None:
        output_accel = min(output_accel, self.standstill_hold_memory)
      if output_accel > self.stop_accel:
        output_accel = min(output_accel, 0.0)
        output_accel -= self.stopping_decel_rate * DT_CTRL
        if soft_hold:
          output_accel = self.stop_accel

      # Once the car is fully stationary, add a little more brake hold without
      # changing the approach-to-stop feel. This prevents the SCC from slowly
      # releasing and re-applying the brakes on a mild downhill or long stop.
      if (standstill_latched or CS.vEgo < 0.05) and not CS.brakePressed:
        # Standstill holding is intentionally independent of StopAccelApply.
        # 래치 중에는 몇 cm 밀려도 계속 건다. 내리막에서 vEgo 가 0.05 를 넘는
        # 순간 홀드가 빠져 더 밀리던 되먹임을 끊는다.
        hold_target = self.standstill_hold_accel
        if output_accel > hold_target:
          output_accel = max(hold_target, output_accel - self.standstill_hold_rate * DT_CTRL)
        self.standstill_hold_memory = output_accel
      self.reset(CS.vEgo)

    elif self.long_control_state == LongCtrlState.starting:
      # Smooth launch: ramp toward StartAccelApply at 2.0 m/s^3 instead of
      # applying the full launch acceleration in one control cycle. At the
      # default control rate this reaches 2.0 m/s^2 in about one second.
      target_start_accel = min(self.start_accel, accel_limits[1])
      self.starting_accel = min(target_start_accel,
                                self.starting_accel + self.starting_ramp_rate * DT_CTRL)
      output_accel = self.starting_accel
      self.reset(CS.vEgo)

    else:
      self.v_pid = v_target_now

      # apilot-c2 low-speed overshoot prevention. Near a planned stop, freeze
      # the integrator so it cannot build an acceleration correction while the
      # car is settling into the final brake hold.
      prevent_overshoot = (not self.CP.stoppingControl and CS.vEgo < 1.5 and
                           v_target_1sec < 0.7 and v_target_1sec < self.v_pid)
      deadzone = interp(CS.vEgo,
                        self.CP.longitudinalTuning.deadzoneBP,
                        self.CP.longitudinalTuning.deadzoneV)
      error = apply_deadzone(self.v_pid - CS.vEgo, deadzone)
      output_accel = self.pid.update(error, speed=CS.vEgo, feedforward=a_target,
                                     freeze_integrator=prevent_overshoot)

      # Keep the starting -> PID handoff bumpless at very low speed. PID may
      # ask for substantially more acceleration than the launch ramp on the
      # first cycle, so only allow the request to rise at the same ramp rate
      # until the car is clearly moving.
      if prev_long_control_state == LongCtrlState.starting and CS.vEgo < 2.0:
        output_accel = min(output_accel,
                           self.last_output_accel + self.starting_ramp_rate * DT_CTRL)
      self.starting_accel = 0.0

      # Back-calculate only a positive saturated integral.  PIDController's
      # regular anti-windup prevents new accumulation, while this removes a
      # correction retained from a previously higher CruiseMax setting.
      raw_control = self.pid.p + self.pid.i + self.pid.d + self.pid.f
      if self.pid.i > 0.0 and raw_control > self.pid.pos_limit:
        self.pid.i = max(0.0, self.pid.pos_limit - self.pid.p - self.pid.d - self.pid.f)
        output_accel = clip(self.pid.p + self.pid.i + self.pid.d + self.pid.f,
                            self.pid.neg_limit, self.pid.pos_limit)
        self.pid.control = output_accel

      if -self.long_coast_band < output_accel < 0.0:
        output_accel = 0.0

    self.last_output_accel = clip(output_accel, accel_limits[0], accel_limits[1])
    return self.last_output_accel, -0.5 if planned_stop else j_target

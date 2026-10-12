from types import SimpleNamespace

from selfdrive.controls.lib.longcontrol import (LEAD_LOST_RELEASE_FRAMES, STANDSTILL_RELEASE_FRAMES,
                                                LongControl, LongCtrlState)


def make_cp():
  return SimpleNamespace(
    enableGasInterceptor=False,
    openpilotLongitudinalControl=True,
    vEgoStopping=0.3,
    vEgoStarting=0.2,
    stopAccel=-0.6,
    stoppingDecelRate=1.2,
    stoppingControl=True,
    longitudinalActuatorDelay=0.5,
    longitudinalActuatorDelayLowerBound=0.5,
    longitudinalActuatorDelayUpperBound=0.5,
    longitudinalTuning=SimpleNamespace(kpBP=[0.], kpV=[1.], kiBP=[0.], kiV=[0.],
                                       kf=1.0, deadzoneBP=[0.], deadzoneV=[0.]),
  )


def make_cs(gas=False, brake=False):
  return SimpleNamespace(vEgo=0.0, brakePressed=brake, gasPressed=gas,
                         cruiseState=SimpleNamespace(standstill=False), buttonEvents=[])


def radar(status=True, d_rel=6.0, v_lead=0.0, is_radar=True):
  return SimpleNamespace(leadOne=SimpleNamespace(status=status, dRel=d_rel, vLead=v_lead, radar=is_radar))


def stopped_controller():
  loc = LongControl(make_cp())
  loc.long_control_state = LongCtrlState.stopping
  return loc


def step(loc, radar_state, v_target_1sec=1.0, cs=None):
  """Feed one 100 Hz frame where the planner asks to launch."""
  return loc._update_standstill_latch(True, cs or make_cs(), 0.0, v_target_1sec, False, radar_state)


def run_until_release(loc, radar_states):
  for i, rs in enumerate(radar_states):
    _, release = step(loc, rs)
    if release:
      return i
  return None


def test_stationary_lead_blocks_gap_filling_launch():
  # 앞차가 서 있는데 플래너가 정지거리를 메우려고 출발을 요구하는 경우.
  loc = stopped_controller()
  assert run_until_release(loc, [radar(d_rel=8.0)] * 1000) is None


def test_radar_lead_departing_releases_after_debounce():
  loc = stopped_controller()
  step(loc, radar(d_rel=6.0))  # latch and record the reference distance
  frame = run_until_release(loc, [radar(d_rel=6.0, v_lead=1.0)] * 200)
  assert frame == STANDSTILL_RELEASE_FRAMES - 1


def test_radar_lead_distance_increase_releases():
  loc = stopped_controller()
  step(loc, radar(d_rel=6.0))
  assert run_until_release(loc, [radar(d_rel=6.5)] * 100) is None
  assert run_until_release(loc, [radar(d_rel=7.2)] * 100) is not None


def test_brief_lead_dropout_does_not_launch():
  loc = stopped_controller()
  step(loc, radar(d_rel=4.0))
  dropout = [radar(status=False)] * (LEAD_LOST_RELEASE_FRAMES - 1)
  assert run_until_release(loc, dropout + [radar(d_rel=4.0)] * 200) is None


def test_lead_gone_long_enough_releases():
  loc = stopped_controller()
  step(loc, radar(d_rel=4.0))
  frame = run_until_release(loc, [radar(status=False)] * 1000)
  assert frame is not None
  assert frame >= LEAD_LOST_RELEASE_FRAMES - 1


def test_vision_lead_noise_needs_larger_departure():
  loc = stopped_controller()
  step(loc, radar(d_rel=5.0, is_radar=False))
  assert run_until_release(loc, [radar(d_rel=6.5, v_lead=0.8, is_radar=False)] * 200) is None
  assert run_until_release(loc, [radar(d_rel=7.5, is_radar=False)] * 200) is not None


def test_no_lead_stop_follows_planner_with_debounce():
  # 신호 대기처럼 앞차 없이 선 경우는 플래너 출발 요구를 0.5 초 뒤 따른다.
  loc = stopped_controller()
  frame = run_until_release(loc, [radar(status=False)] * 200)
  assert frame == STANDSTILL_RELEASE_FRAMES - 1


def test_short_planner_glitch_does_not_launch():
  loc = stopped_controller()
  for _ in range(STANDSTILL_RELEASE_FRAMES - 1):
    _, release = step(loc, None)
    assert not release
  _, release = step(loc, None, v_target_1sec=0.0)
  assert not release
  assert loc.standstill_frames == 0


def test_driver_gas_releases_immediately_even_with_stopped_lead():
  loc = stopped_controller()
  step(loc, radar(d_rel=4.0))
  _, release = step(loc, radar(d_rel=4.0), cs=make_cs(gas=True))
  assert release


def test_lead_tracking_resets_after_release():
  loc = stopped_controller()
  step(loc, radar(d_rel=4.0))
  step(loc, radar(d_rel=4.0), cs=make_cs(gas=True))
  assert not loc.standstill_lead_seen
  assert loc.standstill_lead_ref_d is None

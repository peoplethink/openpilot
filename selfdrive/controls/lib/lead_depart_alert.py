from common.realtime import DT_CTRL
from selfdrive.controls.lib.longcontrol import (LEAD_DEPART_DIST_RADAR, LEAD_DEPART_DIST_VISION,
                                                LEAD_DEPART_SPEED_RADAR, LEAD_DEPART_SPEED_VISION)

# 앞차 출발 알림.
# 내 차가 정차해 있고 바로 앞에 서 있던 차가 출발하면 정차 한 번에 한 번만
# 알린다. 오파 작동 여부와 상관없이 동작한다(신호 대기 중 수동 정차 포함).
EGO_STOPPED_SPEED = 0.1      # m/s 이하면 정차로 본다
LEAD_STOPPED_SPEED = 0.3     # m/s 이하면 앞차도 서 있는 것으로 본다
ARM_TIME = 1.0               # 둘 다 이만큼 서 있어야 알림을 준비한다
CONFIRM_TIME = 0.3           # 출발 판정이 이만큼 이어져야 알린다
MAX_LEAD_DISTANCE = 30.0     # m, 이보다 먼 차는 "바로 앞차"로 보지 않는다


class LeadDepartAlert:
  def __init__(self, dt=DT_CTRL):
    self.arm_frames = max(1, round(ARM_TIME / dt))
    self.confirm_frames = max(1, round(CONFIRM_TIME / dt))
    self.reset()

  def reset(self):
    self.armed = False
    self.fired = False
    self.ref_d = None
    self.stopped_frames = 0
    self.depart_frames = 0

  def update(self, v_ego, gas_pressed, lead):
    """프레임마다 호출. 알림을 띄울 프레임에서만 True 를 돌려준다."""
    if v_ego > EGO_STOPPED_SPEED or gas_pressed:
      # 출발했거나 운전자가 이미 가속 중이면 다음 정차를 위해 초기화한다.
      self.reset()
      return False

    if self.fired:
      return False

    if lead is None or not lead.status or lead.dRel > MAX_LEAD_DISTANCE:
      # 리드가 잠깐 끊긴 것만으로는 알리지 않는다. 준비 상태와 기준거리는
      # 유지해서, 다시 잡혔을 때 멀어져 있으면 그때 알린다.
      self.stopped_frames = 0
      self.depart_frames = 0
      return False

    if not self.armed:
      if lead.vLead < LEAD_STOPPED_SPEED:
        self.stopped_frames += 1
        if self.stopped_frames >= self.arm_frames:
          self.armed = True
          self.ref_d = float(lead.dRel)
      else:
        self.stopped_frames = 0
      return False

    if lead.radar:
      depart_speed, depart_dist = LEAD_DEPART_SPEED_RADAR, LEAD_DEPART_DIST_RADAR
    else:
      depart_speed, depart_dist = LEAD_DEPART_SPEED_VISION, LEAD_DEPART_DIST_VISION
    departing = lead.vLead > depart_speed or lead.dRel - self.ref_d > depart_dist
    self.depart_frames = self.depart_frames + 1 if departing else 0

    if self.depart_frames >= self.confirm_frames:
      self.fired = True
      return True
    return False

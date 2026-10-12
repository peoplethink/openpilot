from types import SimpleNamespace

from selfdrive.controls.lib.lead_depart_alert import LeadDepartAlert


def lead(d_rel=6.0, v_lead=0.0, status=True, is_radar=True):
  return SimpleNamespace(status=status, dRel=d_rel, vLead=v_lead, radar=is_radar)


def run(alert, frames, v_ego=0.0, gas=False):
  """Return the indexes of frames where the alert fired."""
  return [i for i, ld in enumerate(frames) if alert.update(v_ego, gas, ld)]


def stopped_behind_lead(seconds=1.5, d_rel=6.0, is_radar=True):
  return [lead(d_rel=d_rel, is_radar=is_radar)] * int(seconds * 100)


def test_fires_once_when_stopped_lead_pulls_away():
  alert = LeadDepartAlert()
  frames = stopped_behind_lead() + [lead(d_rel=6.0, v_lead=1.0)] * 200
  fired = run(alert, frames)
  assert len(fired) == 1
  assert fired[0] == 150 + 29  # 0.3 s confirmation after the lead starts moving


def test_distance_increase_alone_fires():
  alert = LeadDepartAlert()
  fired = run(alert, stopped_behind_lead() + [lead(d_rel=7.3)] * 50)
  assert len(fired) == 1


def test_no_alert_while_lead_stays_stopped():
  alert = LeadDepartAlert()
  assert run(alert, stopped_behind_lead(seconds=20)) == []


def test_short_stop_does_not_arm():
  # 1 초도 안 서 있었으면(서행 정체 중 잠깐 멈춤) 알리지 않는다.
  alert = LeadDepartAlert()
  assert run(alert, stopped_behind_lead(seconds=0.5) + [lead(v_lead=1.0)] * 100) == []


def test_lead_dropout_alone_does_not_fire():
  alert = LeadDepartAlert()
  assert run(alert, stopped_behind_lead() + [lead(status=False)] * 500) == []


def test_lead_reacquired_farther_away_fires():
  alert = LeadDepartAlert()
  frames = stopped_behind_lead() + [lead(status=False)] * 100 + [lead(d_rel=9.0, v_lead=0.2)] * 50
  assert len(run(alert, frames)) == 1


def test_vision_noise_needs_bigger_change():
  alert = LeadDepartAlert()
  frames = stopped_behind_lead(is_radar=False) + [lead(d_rel=7.5, v_lead=0.8, is_radar=False)] * 100
  assert run(alert, frames) == []


def test_no_alert_when_driver_is_already_going():
  alert = LeadDepartAlert()
  assert run(alert, stopped_behind_lead()) == []
  assert run(alert, [lead(v_lead=1.0)] * 100, gas=True) == []


def test_rearms_for_the_next_stop():
  alert = LeadDepartAlert()
  assert len(run(alert, stopped_behind_lead() + [lead(v_lead=1.0)] * 50)) == 1
  run(alert, [lead(v_lead=1.0)], v_ego=3.0)  # drive off
  assert len(run(alert, stopped_behind_lead() + [lead(v_lead=1.0)] * 50)) == 1


def test_far_lead_is_ignored():
  alert = LeadDepartAlert()
  assert run(alert, stopped_behind_lead(d_rel=40.0) + [lead(d_rel=40.0, v_lead=1.0)] * 100) == []

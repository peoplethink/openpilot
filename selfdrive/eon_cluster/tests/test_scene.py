from types import SimpleNamespace

from selfdrive.eon_cluster.scene import camera_lane_position, final_lateral_path, reconcile_lane_position


def polyline(xs, ys):
  return SimpleNamespace(x=xs, y=ys)


def test_final_lateral_path_uses_optimized_mpc_xy_and_model_z():
  lateral = SimpleNamespace(
    mpcSolutionValid=True,
    dPathPoints=[0.0, 9.0, 18.0],
    mpcPathX=[0.0, 0.9, 3.7],
    mpcPathY=[0.0, 0.4, 1.8],
  )
  model = SimpleNamespace(
    position=SimpleNamespace(z=[0.0, 0.1, 0.2]),
  )
  assert final_lateral_path(lateral, model, [0.0, 0.1, 0.4]) == [
    [0.0, 0.0, 0.0], [0.9, 0.4, 0.1], [3.7, 1.8, 0.2],
  ]


def test_final_lateral_path_rejects_invalid_or_missing_optimized_mpc_path():
  model = SimpleNamespace(position=SimpleNamespace(z=[]))
  invalid = SimpleNamespace(mpcSolutionValid=False,
                            mpcPathX=[0.0, 1.0], mpcPathY=[0.0, 1.0])
  missing_solution = SimpleNamespace(mpcSolutionValid=True,
                                     dPathPoints=[0.0, 1.0])
  assert final_lateral_path(invalid, model, [0.0, 1.0]) == []
  assert final_lateral_path(missing_solution, model, [0.0, 1.0]) == []


def lane_position_model(left_edge, right_edge, lane_probs=(0.8, 0.9, 0.9, 0.8),
                        edge_stds=(0.1, 0.1)):
  xs = [0.0, 10.0, 20.0, 30.0]
  return SimpleNamespace(
    laneLines=[polyline(xs, [5.4] * 4), polyline(xs, [1.8] * 4),
               polyline(xs, [-1.8] * 4), polyline(xs, [-5.4] * 4)],
    laneLineProbs=list(lane_probs),
    roadEdges=[polyline(xs, [left_edge] * 4), polyline(xs, [right_edge] * 4)],
    roadEdgeStds=list(edge_stds),
  )


def test_camera_lane_position_tracks_left_middle_and_right_lane():
  assert camera_lane_position(lane_position_model(2.5, -9.0))["cur"] == 1
  middle = camera_lane_position(lane_position_model(5.4, -5.4))
  assert middle["n"] == 3 and middle["cur"] == 2
  assert camera_lane_position(lane_position_model(9.0, -2.5))["cur"] == 3


def test_camera_lane_position_accepts_edge_order_and_rejects_weak_geometry():
  reversed_edges = lane_position_model(-5.4, 5.4)
  assert camera_lane_position(reversed_edges) == {
    "n": 3, "cur": 2, "confidence": 0.9, "laneWidth": 3.6,
    "leftFrac": 0.0, "rightFrac": 0.0,
  }
  assert camera_lane_position(lane_position_model(
    5.4, -5.4, lane_probs=(0.8, 0.2, 0.9, 0.8))) is None
  assert camera_lane_position(lane_position_model(
    5.4, -5.4, edge_stds=(0.8, 0.1))) is None


def test_camera_lane_position_reports_lane_fractions():
  # road edge 1.4 lanes left of the ego lane, 0.0 lanes right
  position = camera_lane_position(lane_position_model(1.8 + 3.6 * 1.4, -1.8))
  assert position["n"] == 2 and position["cur"] == 2
  assert abs(position["leftFrac"] - 0.4) < 0.01
  assert position["rightFrac"] == 0.0


def test_reconcile_keeps_matching_count():
  position = {"n": 3, "cur": 2, "leftFrac": 0.0, "rightFrac": 0.0}
  assert reconcile_lane_position(position, 3) == {**position, "source": "route"}


def test_reconcile_drops_right_shoulder_counted_as_a_lane():
  # 오른쪽 갓길이 0.6 차로폭 → 반올림으로 차로 하나가 생긴 경우.
  position = {"n": 3, "cur": 2, "leftFrac": 0.0, "rightFrac": 0.6}
  reconciled = reconcile_lane_position(position, 2)
  assert reconciled["n"] == 2 and reconciled["cur"] == 2


def test_reconcile_drops_left_median_counted_as_a_lane():
  position = {"n": 3, "cur": 2, "leftFrac": 0.7, "rightFrac": 0.0}
  reconciled = reconcile_lane_position(position, 2)
  assert reconciled["n"] == 2 and reconciled["cur"] == 1


def test_reconcile_fails_closed_when_ambiguous_or_far_off():
  both = {"n": 3, "cur": 2, "leftFrac": 0.6, "rightFrac": 0.6}
  whole_lane = {"n": 3, "cur": 2, "leftFrac": 0.95, "rightFrac": 0.0}
  assert reconcile_lane_position(both, 2) is None
  assert reconcile_lane_position(whole_lane, 2) is None
  assert reconcile_lane_position({"n": 5, "cur": 2}, 2) is None
  assert reconcile_lane_position(None, 2) is None
  assert reconcile_lane_position({"n": 2, "cur": 1}, 0) is None
  assert reconcile_lane_position({"n": 2, "cur": 1}, "bad") is None

import math
import statistics


def _field(message, name, default=None):
  try:
    return getattr(message, name)
  except Exception:
    return default


def _finite_float(value, default=0.0):
  try:
    number = float(value)
  except (TypeError, ValueError):
    return default
  return number if math.isfinite(number) else default


def _near_line_y(line, maximum_x=30.0):
  """Robust near-field lateral position for HUD-only lane placement."""
  xs = list(_field(line, "x", []) or [])
  ys = list(_field(line, "y", []) or [])
  values = [_finite_float(y, float("nan")) for x, y in zip(xs, ys)
            if math.isfinite(_finite_float(x, float("nan"))) and
            0.0 <= _finite_float(x) <= maximum_x and
            math.isfinite(_finite_float(y, float("nan")))]
  return statistics.median(values) if len(values) >= 2 else None


ROAD_EDGE_STD_MAX = 0.5
# 갓길·중앙분리대 유령차로 판정 여유. navigation_noo 의
# NavigationLaneChangeController.PHANTOM_LANE_MIN_ERROR 와 같은 값.
PHANTOM_LANE_MIN_ERROR = 0.10


def camera_lane_position(model):
  """Estimate the ego lane from modelV2 lane lines and road edges.

  modelV2 does not publish a semantic "lane 2" value. It does publish the two
  ego-lane boundaries and both road edges, so the number of lane-width spaces
  between them gives a conservative camera-relative lane index. This metadata
  is display-only and never feeds lateral control.
  """
  lanes = list(_field(model, "laneLines", []) or [])
  edges = list(_field(model, "roadEdges", []) or [])
  lane_probs = list(_field(model, "laneLineProbs", []) or [])
  edge_stds = list(_field(model, "roadEdgeStds", []) or [])
  if len(lanes) < 3 or len(edges) < 2:
    return None

  try:
    inner = [_near_line_y(lanes[1]), _near_line_y(lanes[2])]
    road = [_near_line_y(edges[0]), _near_line_y(edges[1])]
    inner_conf = min(_finite_float(lane_probs[1], 0.0), _finite_float(lane_probs[2], 0.0))
    # roadEdgeStds is a standard deviation in metres, not a probability, so it
    # is compared against a distance threshold (same rule as NOO control).
    edge_std = max(_finite_float(edge_stds[0], 9.9), _finite_float(edge_stds[1], 9.9))
  except (IndexError, TypeError):
    return None
  if any(value is None for value in inner + road) or inner_conf < 0.45 or \
     edge_std > ROAD_EDGE_STD_MAX:
    return None

  lane_left, lane_right = max(inner), min(inner)
  road_left, road_right = max(road), min(road)
  lane_width = lane_left - lane_right
  if not 2.5 <= lane_width <= 4.5 or road_left < lane_left or road_right > lane_right:
    return None

  left_ratio = max(0.0, road_left - lane_left) / lane_width
  right_ratio = max(0.0, lane_right - road_right) / lane_width
  left_lanes = int(round(left_ratio))
  right_lanes = int(round(right_ratio))
  total = 1 + left_lanes + right_lanes
  current = 1 + left_lanes
  if not 1 <= total <= 8 or not 1 <= current <= total:
    return None
  return {
    "n": total,
    "cur": current,
    "confidence": round(inner_conf, 2),
    "laneWidth": round(lane_width, 2),
    # 그 쪽 여유폭이 차로폭의 몇 배인지의 소수부. reconcile_lane_position 이
    # 갓길·중앙분리대가 차로 하나로 잡혔는지 판정할 때 쓴다(NOO 와 같은 기준).
    "leftFrac": round(left_ratio - math.floor(left_ratio), 3),
    "rightFrac": round(right_ratio - math.floor(right_ratio), 3),
  }


def reconcile_lane_position(lane_position, route_lane_count):
  """카메라 차로 추정을 TMAP 차로 수에 맞춘다(HUD 표시 전용).

  한국 도로는 갓길이 있어 카메라가 차로를 하나 더 세는 일이 흔하다.
  NOO 차로 판정(NavigationLaneChangeController._resolve_current_lane)과
  같은 규칙으로, 카메라가 정확히 하나 더 셌고 어느 쪽이 유령차로인지
  분명할 때만 그 쪽을 빼서 현재 차로를 고친다. 판정할 수 없으면 None 을
  돌려주고, 호출하는 쪽은 카메라 추정을 그대로 쓴다.
  """
  if not isinstance(lane_position, dict):
    return None
  try:
    route_count = int(route_lane_count)
    camera_count = int(lane_position.get("n", 0))
    current = int(lane_position.get("cur", 0))
  except (TypeError, ValueError):
    return None
  if not 1 <= route_count <= 8 or not 1 <= current <= camera_count:
    return None

  if camera_count == route_count:
    resolved = current
  elif camera_count - route_count == 1:
    def phantom(frac):
      return 0.5 <= frac <= 1.0 - PHANTOM_LANE_MIN_ERROR

    phantom_left = phantom(_finite_float(lane_position.get("leftFrac"), 0.0))
    phantom_right = phantom(_finite_float(lane_position.get("rightFrac"), 0.0))
    if phantom_left == phantom_right:
      return None
    resolved = current - 1 if phantom_left else current
  else:
    return None

  if not 1 <= resolved <= route_count:
    return None
  reconciled = dict(lane_position)
  reconciled["n"] = route_count
  reconciled["cur"] = resolved
  reconciled["source"] = "route"
  return reconciled


def final_lateral_path(lateral_plan, model, time_indices, limit=33):
  """Build the optimized MPC path and retain model road elevation.

  dPathPoints is only the reference passed into the solver. LateralPlanner
  publishes the solver's x/y state trajectory as mpcPathX/mpcPathY. Keep model
  Z at the same shooting-node index so the renderer retains road elevation.
  """
  if not bool(_field(lateral_plan, "mpcSolutionValid", False)):
    return []
  raw_xs = _field(lateral_plan, "mpcPathX", [])
  raw_ys = _field(lateral_plan, "mpcPathY", [])
  xs = list(raw_xs) if raw_xs is not None else []
  ys = list(raw_ys) if raw_ys is not None else []
  if not xs or not ys:
    return []

  position = _field(model, "position")
  zs = list(_field(position, "z", []) or [])
  count = min(len(xs), len(ys), len(time_indices), int(limit))
  if count < 2:
    return []
  return [[round(_finite_float(xs[i]), 2),
           round(_finite_float(ys[i]), 2),
           round(_finite_float(zs[i] if i < len(zs) else 0.0), 2)]
          for i in range(count)]

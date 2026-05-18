"""Helpers for the optional ``splitLines`` and ``visualSplitGroupId``
fields carried by items of ``integration.stack``.
"""

from __future__ import annotations

import math
import re
from typing import Dict, Iterable, List, Optional


_BOUND_EPSILON = 1e-9
_DEDUP_EPSILON = 1e-6

# Mirrors the frontend validation regex /[^A-Za-z0-9._-]/g.
_GROUP_ID_RE = re.compile(r'^[A-Za-z0-9._-]+$')

def _coerce_float(value) -> Optional[float]:
    if value is None or isinstance(value, bool):
        return None
    try:
        as_float = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(as_float) or math.isinf(as_float):
        return None
    return as_float


def _format_number(value: float) -> str:
    text = repr(float(value))
    if text.endswith('.0'):
        return text[:-2]
    return text


def sanitize_split_lines(
    values,
    x_lower: Optional[float] = None,
    x_upper: Optional[float] = None,
) -> Optional[List[float]]:
    if values is None:
        return None
    if isinstance(values, (str, bytes)):
        return None
    if not isinstance(values, Iterable):
        return None

    lo: Optional[float] = None
    hi: Optional[float] = None
    if x_lower is not None and x_upper is not None:
        lo, hi = sorted((float(x_lower), float(x_upper)))

    cleaned: List[float] = []
    for raw in values:
        as_float = _coerce_float(raw)
        if as_float is None:
            continue
        if lo is not None and (as_float < lo - _BOUND_EPSILON or as_float > hi + _BOUND_EPSILON):
            continue
        cleaned.append(as_float)

    if not cleaned:
        return None

    cleaned.sort()
    deduped: List[float] = [cleaned[0]]
    for value in cleaned[1:]:
        if abs(value - deduped[-1]) > _DEDUP_EPSILON:
            deduped.append(value)

    return deduped or None


def format_split_lines(values) -> Optional[str]:
    """Render ``splitLines`` as the JCAMP token ``"v1;v2;v3"`` or ``None``."""
    if not values:
        return None
    return ';'.join(_format_number(v) for v in values)


def parse_split_lines(text) -> Optional[List[float]]:
    """Inverse of :func:`format_split_lines`. Tolerates surrounding quotes."""
    if text is None:
        return None
    if not isinstance(text, str):
        return None
    cleaned = text.strip().strip('"').strip("'").strip()
    if not cleaned:
        return None
    parts = [chunk.strip() for chunk in cleaned.split(';') if chunk.strip()]
    out: List[float] = []
    for part in parts:
        as_float = _coerce_float(part)
        if as_float is not None:
            out.append(as_float)
    return out or None


def sanitize_group_id(value) -> Optional[str]:
    """Return a cleaned ``visualSplitGroupId`` or ``None``."""
    if value is None:
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, str):
        cleaned = value.strip()
    else:
        try:
            cleaned = str(value).strip()
        except Exception:
            return None
    if not cleaned:
        return None
    if not _GROUP_ID_RE.match(cleaned):
        return None
    return cleaned


def format_observed_integrals_groups_table(items: Iterable[dict]) -> List[str]:
    out: List[str] = []
    for idx, item in enumerate(items):
        if not isinstance(item, dict):
            continue
        group_id = sanitize_group_id(item.get('visualSplitGroupId'))
        if group_id is None:
            continue
        out.append('{}, {}\n'.format(idx, group_id))
    return out


def parse_observed_integrals_groups_table(text: str) -> Dict[int, str]:
    if not text:
        return {}
    mapping: Dict[int, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip().lstrip('(').rstrip(')')
        if not line:
            continue
        parts = [c.strip() for c in line.split(',')]
        if len(parts) < 2:
            continue
        try:
            idx = int(parts[0])
        except (TypeError, ValueError):
            continue
        group_id_raw = ','.join(parts[1:]).strip().strip('"').strip("'")
        group_id = sanitize_group_id(group_id_raw)
        if group_id is None:
            continue
        mapping[idx] = group_id
    return mapping


def parse_observed_integrals_table(
    text: str,
    groups_text: Optional[str] = None,
) -> List[dict]:
    if not text:
        return []

    rows: List[dict] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        line = line.lstrip('(').rstrip(')')
        cols = [c.strip() for c in line.split(',')]
        if len(cols) < 3:
            continue
        xL = _coerce_float(cols[0])
        xU = _coerce_float(cols[1])
        area = _coerce_float(cols[2])
        if xL is None or xU is None or area is None:
            continue
        item: dict = {'xL': xL, 'xU': xU, 'area': area}
        if len(cols) >= 4:
            abs_area = _coerce_float(cols[3])
            if abs_area is not None:
                item['absoluteArea'] = abs_area
        if len(cols) >= 5:
            split_lines = parse_split_lines(cols[4])
            if split_lines:
                item['splitLines'] = split_lines
        rows.append(item)

    if groups_text:
        groups = parse_observed_integrals_groups_table(groups_text)
        for idx, group_id in groups.items():
            if 0 <= idx < len(rows):
                rows[idx]['visualSplitGroupId'] = group_id
    return rows

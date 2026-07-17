def is_valid_integration_item(itg):
    if not isinstance(itg, dict):
        return False
    if itg.get('xL') is None or itg.get('xU') is None:
        return False
    return itg.get('area') is not None


def filter_valid_integrations(stack):
    return [itg for itg in (stack or []) if is_valid_integration_item(itg)]


def integration_uses_auc_column(items, is_hplc_uv_vis=False):
    if is_hplc_uv_vis:
        return True
    return any('absoluteArea' in itg for itg in items)


def format_integral_line(itg, ref_area_factor, ref_shift, use_auc_column):
    absolute_area = itg.get('absoluteArea') or 0
    area_stored = float(itg.get('area') or 0) * ref_area_factor
    x_l = itg['xL'] - ref_shift
    x_u = itg['xU'] - ref_shift
    if use_auc_column:
        return '({}, {}, {}, {})\n'.format(x_l, x_u, area_stored, absolute_area)
    return '({}, {}, {})\n'.format(x_l, x_u, area_stored)


def build_integration_lines(items, ref_area_factor, ref_shift, use_auc_column):
    return [
        format_integral_line(itg, ref_area_factor, ref_shift, use_auc_column)
        for itg in items
        if is_valid_integration_item(itg)
    ]


def serialize_integration_stack(itg_stack, ref_area_factor, ref_shift, is_hplc_uv_vis=False):
    if not itg_stack:
        return []
    if isinstance(itg_stack[0], str):
        return itg_stack
    use_auc_column = integration_uses_auc_column(itg_stack, is_hplc_uv_vis)
    return build_integration_lines(itg_stack, ref_area_factor, ref_shift, use_auc_column)


def valid_visual_split_indices(items):
    valid_indices = set()
    idx = 0
    while idx < len(items):
        group_id = items[idx].get('visualSplitGroupId')
        if not group_id:
            idx += 1
            continue
        end = idx + 1
        while end < len(items) and items[end].get('visualSplitGroupId') == group_id:
            end += 1
        if end - idx >= 2:
            valid_indices.update(range(idx, end))
        idx = end if end > idx + 1 else idx + 1
    return valid_indices


def group_visual_splits(items, ref_shift=0):
    groups = []
    idx = 0
    while idx < len(items):
        group_id = items[idx].get('visualSplitGroupId')
        if not group_id:
            groups.append({'items': items[idx:idx + 1], 'split_xs': []})
            idx += 1
            continue
        end = idx + 1
        while end < len(items) and items[end].get('visualSplitGroupId') == group_id:
            end += 1
        if end - idx >= 2:
            split_xs = [
                items[i]['xU'] - ref_shift
                for i in range(idx, end - 1)
            ]
            groups.append({'items': items[idx:end], 'split_xs': split_xs})
        else:
            groups.append({'items': [items[idx]], 'split_xs': []})
        idx = end if end > idx + 1 else idx + 1
    return groups


def build_integration_groups_lines(items):
    valid_indices = valid_visual_split_indices(items)
    if not valid_indices:
        return []
    table = []
    for idx, itg in enumerate(items):
        if idx not in valid_indices:
            continue
        group_id = itg.get('visualSplitGroupId')
        if group_id:
            table.append('{}, {}\n'.format(idx, group_id))
    return table

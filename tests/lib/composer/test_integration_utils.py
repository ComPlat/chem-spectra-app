from chem_spectra.lib.composer.base import is_metadata_to_be_ignored
from chem_spectra.lib.composer.integration_utils import (
    build_integration_groups_lines,
    build_integration_lines,
    filter_valid_integrations,
    group_visual_splits,
    integration_uses_auc_column,
    is_valid_integration_item,
    serialize_integration_stack,
)


def test_is_valid_integration_item():
    assert is_valid_integration_item({'xL': 0, 'xU': 1, 'area': 1.0})
    assert not is_valid_integration_item({'xL': 0, 'xU': 1, 'area': None})
    assert not is_valid_integration_item('not a dict')


def test_filter_valid_integrations():
    stack = [
        {'xL': 0, 'xU': 1, 'area': 1.0},
        {'xL': 2, 'xU': 3, 'area': None},
    ]
    assert len(filter_valid_integrations(stack)) == 1


def test_integration_uses_auc_column():
    items = [{'xL': 0, 'xU': 1, 'area': 1.0}]
    assert integration_uses_auc_column(items, is_hplc_uv_vis=True)
    assert not integration_uses_auc_column(items, is_hplc_uv_vis=False)
    assert integration_uses_auc_column(
        [{'xL': 0, 'xU': 1, 'area': 1.0, 'absoluteArea': 10}],
        is_hplc_uv_vis=False,
    )


def test_build_integration_lines_nmr_three_columns():
    items = [{'xL': 0.0, 'xU': 4.0, 'area': 4.0}]
    lines = build_integration_lines(items, ref_area_factor=0.5, ref_shift=0, use_auc_column=False)
    assert lines == ['(0.0, 4.0, 2.0)\n']


def test_build_integration_lines_hplc_four_columns():
    items = [{'xL': 0.0, 'xU': 4.0, 'area': 4.0, 'absoluteArea': 40.0}]
    lines = build_integration_lines(items, ref_area_factor=1.0, ref_shift=0, use_auc_column=True)
    assert lines == ['(0.0, 4.0, 4.0, 40.0)\n']


def test_group_visual_splits():
    items = [
        {'xL': 1.0, 'xU': 2.0, 'visualSplitGroupId': 'g1'},
        {'xL': 2.0, 'xU': 3.0, 'visualSplitGroupId': 'g1'},
        {'xL': 5.0, 'xU': 6.0},
    ]
    groups = group_visual_splits(items, ref_shift=0)
    assert len(groups) == 2
    assert groups[0]['split_xs'] == [2.0]


def test_build_integration_groups_lines():
    items = [
        {'xL': 0, 'xU': 4, 'area': 1, 'visualSplitGroupId': 'g1'},
        {'xL': 4, 'xU': 10, 'area': 2, 'visualSplitGroupId': 'g1'},
    ]
    assert build_integration_groups_lines(items) == ['0, g1\n', '1, g1\n']


def test_serialize_integration_stack_from_strings():
    assert serialize_integration_stack(['(1, 2, 3)\n'], 1, 0) == ['(1, 2, 3)\n']


def test_original_metadata_ignores_chemspectra_keys():
    assert is_metadata_to_be_ignored('$OBSERVEDINTEGRALS')
    assert is_metadata_to_be_ignored('$OBSERVEDINTEGRALSGROUPS')
    assert is_metadata_to_be_ignored('$CSITAREA')

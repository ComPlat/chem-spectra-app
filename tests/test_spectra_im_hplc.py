from test_spectra_peaks import (
    __generated_composer,
    __generated_peaks_meta,
    __target_peaks_meta,
)

HPLC_jdx = 'hplc/chromatogram.jdx'

empty_multiplicity = '{"stack":[],"shift":0,"smExtext":false}'

hplc_visual_split_params = {
    'integration': (
        '{"stack":['
        '{"xL":0.0,"xU":4.0,"area":4.0,"absoluteArea":40.0,"visualSplitGroupId":"vsg-abc123-1"},'
        '{"xL":4.0,"xU":10.0,"area":6.0,"absoluteArea":60.0,"visualSplitGroupId":"vsg-abc123-1"}'
        '],"refArea":1.0,"refFactor":1.0,"shift":0,"edited":true}'
    ),
    'multiplicity': empty_multiplicity,
}

hplc_simple_split_params = {
    'integration': (
        '{"stack":['
        '{"xL":0.0,"xU":4.0,"area":2.0,"absoluteArea":20.0},'
        '{"xL":4.0,"xU":10.0,"area":3.0,"absoluteArea":30.0}'
        '],"refArea":1.0,"refFactor":1.0,"shift":0,"edited":true}'
    ),
    'multiplicity': empty_multiplicity,
}

hplc_csit_params = {
    'integration': (
        '{"stack":[{"xL":0.0,"xU":4.0,"area":2.0,"absoluteArea":20.0}],'
        '"refArea":1.5,"refFactor":2,"shift":0,"edited":true}'
    ),
    'multiplicity': empty_multiplicity,
}

hplc_empty_params = {
    'integration': (
        '{"stack":[],"refArea":1.0,"refFactor":1.0,"shift":0,'
        '"edited":true,"originStack":true}'
    ),
    'multiplicity': empty_multiplicity,
}

hplc_orphan_split_params = {
    'integration': (
        '{"stack":['
        '{"xL":0.0,"xU":4.0,"area":2.0,"visualSplitGroupId":"vsg-orphan-1"},'
        '{"xL":4.0,"xU":10.0,"area":3.0,"absoluteArea":30.0}'
        '],"refArea":1.0,"refFactor":1.0,"shift":0,"edited":true}'
    ),
    'multiplicity': empty_multiplicity,
}

hplc_mixed_groups_params = {
    'integration': (
        '{"stack":['
        '{"xL":0.0,"xU":2.0,"area":1.0,"absoluteArea":10.0},'
        '{"xL":2.0,"xU":4.0,"area":2.0,"absoluteArea":20.0,"visualSplitGroupId":"g1"},'
        '{"xL":4.0,"xU":6.0,"area":3.0,"absoluteArea":30.0,"visualSplitGroupId":"g1"}'
        '],"refArea":1.0,"refFactor":1.0,"shift":0,"edited":true}'
    ),
    'multiplicity': empty_multiplicity,
}

hplc_null_area_params = {
    'integration': (
        '{"stack":[{"xL":1.0,"xU":2.0,"area":null}],'
        '"refArea":1.0,"refFactor":1.0,"shift":0,"edited":true}'
    ),
    'multiplicity': empty_multiplicity,
}


def test_im_params_meta_hplc_visual_split():
    meta_content = __generated_peaks_meta(HPLC_jdx, hplc_visual_split_params)
    meta_target = __target_peaks_meta('im/im_meta_hplc_visual_split')
    assert meta_content == meta_target


def test_im_params_meta_hplc_simple_split():
    meta_content = __generated_peaks_meta(HPLC_jdx, hplc_simple_split_params)
    meta_target = __target_peaks_meta('im/im_meta_hplc_simple_split')
    assert meta_content == meta_target


def test_im_params_meta_hplc_csit():
    meta_content = __generated_peaks_meta(HPLC_jdx, hplc_csit_params)
    meta_target = __target_peaks_meta('im/im_meta_hplc_csit')
    assert meta_content == meta_target


def test_im_params_meta_hplc_empty_integrations():
    meta_content = __generated_peaks_meta(HPLC_jdx, hplc_empty_params)
    meta_target = __target_peaks_meta('im/im_meta_hplc_empty')
    assert meta_content == meta_target


def test_im_params_meta_hplc_orphan_visual_split():
    meta_content = __generated_peaks_meta(HPLC_jdx, hplc_orphan_split_params)
    meta_target = __target_peaks_meta('im/im_meta_hplc_orphan_split')
    assert meta_content == meta_target


def test_im_params_meta_hplc_mixed_visual_split_groups():
    meta_content = __generated_peaks_meta(HPLC_jdx, hplc_mixed_groups_params)
    meta_target = __target_peaks_meta('im/im_meta_hplc_mixed_groups')
    assert meta_content == meta_target


def test_hplc_no_multiplet_headers():
    composer = __generated_composer(HPLC_jdx, hplc_visual_split_params)
    joined = composer._join_meta(composer.meta)
    assert '##$OBSERVEDMULTIPLETS' not in joined


def test_hplc_null_area_clears_integrations():
    composer = __generated_composer(HPLC_jdx, hplc_null_area_params)
    assert composer.gen_integration_info() == []


def test_hplc_tf_jcamp_meta_contains_only_strings():
    composer = __generated_composer(HPLC_jdx, hplc_visual_split_params)
    joined = composer._join_meta(composer.meta)
    assert isinstance(joined, str)
    assert '##$OBSERVEDINTEGRALS' in joined

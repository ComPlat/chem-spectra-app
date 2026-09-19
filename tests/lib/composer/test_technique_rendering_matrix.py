"""Per-technique rendering contract.

For every key in data_type.json this asserts the four things the sixteen
`is_*` booleans currently decide between them: peak threshold, x-axis
orientation, x-label style and y-label style.

The axis assertions capture what the composer actually draws, by spying on
`plt.xlim`/`plt.xlabel`/`plt.ylabel`, rather than re-deriving the condition
in `NIComposer.tf_img`. Re-deriving it would make the test a copy of the
code and it would agree with any regression.

Expected values were extracted by execution, and DSC and LC/MS carry the
corrected forward orientation this branch fixes.

The unrecognised-datatype case is included: since #291 such a file takes
the generic curve path instead of raising from __index_target.
"""

import json
import os

import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
# NIComposer selects the Agg backend on import, so import pyplot after it
from chem_spectra.lib.composer.ni import NIComposer
import chem_spectra.lib.composer.ni as ni_module

SOURCE = './tests/fixtures/source/hplc/chromatogram.jdx'
HEADER = '##DATA TYPE=HPLC UV/VIS SPECTRUM'

# key -> (threshold, orientation, x-label style, y-label style)
EXPECTED = {
    'NMR':                               (0.005, 'reversed', 'chemical_shift', 'intensity'),
    'INFRARED':                          (0.93,  'reversed', 'generic',        'generic'),
    'RAMAN':                             (0.07,  'reversed', 'generic',        'generic'),
    'MS':                                (0.05,  'reversed', 'generic',        'generic'),
    'HPLC UVVIS':                        (0.05,  'forward',  'generic',        'generic'),
    'UVVIS':                             (0.05,  'forward',  'generic',        'generic'),
    'THERMOGRAVIMETRIC ANALYSIS':        (1.05,  'forward',  'generic',        'generic'),
    'X-RAY DIFFRACTION':                 (1.00,  'forward',  'xrd',            'generic'),
    'CYCLIC VOLTAMMETRY':                (1.00,  'forward',  'raw',            'raw'),
    'SIZE EXCLUSION CHROMATOGRAPHY':     (0.5,   'forward',  'generic',        'generic'),
    'CIRCULAR DICHROISM SPECTROSCOPY':   (1.00,  'forward',  'generic',        'generic'),
    'SORPTION-DESORPTION MEASUREMENT':   (1.00,  'forward',  'generic',        'generic'),
    'Emissions':                         (0.5,   'forward',  'generic',        'generic'),
    'DLS ACF':                           (1.05,  'forward',  'generic',        'generic'),
    'DLS intensity':                     (1.00,  'forward',  'generic',        'generic'),
    # forward: was reversed because is_dsc was never
    # added to the orientation chain. TGA, the same family, was forward.
    'DIFFERENTIAL SCANNING CALORIMETRY': (1.05,  'forward',  'generic',        'generic'),
    'GAS CHROMATOGRAPHY':                (0.5,   'forward',  'generic',        'generic'),
    # forward: a chromatogram runs forward in
    # retention time. It had no flag and inherited the unknown default.
    'LC/MS':                             (0.5,   'forward',  'generic',        'generic'),
}


def _mapping_keys():
    directory = os.path.dirname(
        __import__('chem_spectra.lib.converter.jcamp.base', fromlist=['base']).__file__
    )
    with open(os.path.join(directory, 'data_type.json')) as handle:
        return json.load(handle)['datatypes']


@pytest.fixture
def render():
    """Render a spectrum of the given datatype; report what it drew."""
    template = open(SOURCE).read()
    saved = (ni_module.plt.xlim, ni_module.plt.xlabel, ni_module.plt.ylabel)

    def _render(datatype, tmp_path):
        target = tmp_path / 'probe.jdx'
        target.write_text(template.replace(HEADER, '##DATA TYPE=' + datatype, 1))

        drawn = {}
        real_xlim, real_xlabel, real_ylabel = saved

        def spy_xlim(*args, **kwargs):
            if args:
                drawn['xlim'] = args
            return real_xlim(*args, **kwargs)

        def spy_xlabel(text, **kwargs):
            drawn['xlabel'] = text
            return real_xlabel(text, **kwargs)

        def spy_ylabel(text, **kwargs):
            drawn['ylabel'] = text
            return real_ylabel(text, **kwargs)

        ni_module.plt.xlim = spy_xlim
        ni_module.plt.xlabel = spy_xlabel
        ni_module.plt.ylabel = spy_ylabel
        try:
            converter = JcampNIConverter(JcampBaseConverter(str(target)))
            composer = NIComposer(converter)
            composer.tf_img().close()
        finally:
            ni_module.plt.xlim, ni_module.plt.xlabel, ni_module.plt.ylabel = saved

        low, high = drawn['xlim']
        drawn['orientation'] = 'forward' if low < high else 'reversed'
        drawn['threshold'] = converter.threshold
        return drawn

    yield _render
    ni_module.plt.xlim, ni_module.plt.xlabel, ni_module.plt.ylabel = saved


def _x_style(label):
    if label.startswith('Chemical shift'):
        return 'chemical_shift'
    if ', WL=' in label:
        return 'xrd'
    if label.startswith('X ('):
        return 'generic'
    return 'raw'


def _y_style(label):
    if label.startswith('Intensity'):
        return 'intensity'
    if label.startswith('Y ('):
        return 'generic'
    return 'raw'


@pytest.mark.parametrize('key', sorted(EXPECTED))
def test_technique_rendering_contract(key, render, tmp_path):
    datatype = _mapping_keys()[key][0]
    drawn = render(datatype, tmp_path)
    threshold, orientation, x_style, y_style = EXPECTED[key]

    assert drawn['threshold'] == threshold, 'threshold for %s' % key
    assert drawn['orientation'] == orientation, 'x-orientation for %s' % key
    assert _x_style(drawn['xlabel']) == x_style, 'x-label style for %s (%r)' % (key, drawn['xlabel'])
    assert _y_style(drawn['ylabel']) == y_style, 'y-label style for %s (%r)' % (key, drawn['ylabel'])


def test_matrix_covers_every_mapped_technique():
    """If a key is added to data_type.json, this matrix must gain a row."""
    assert set(_mapping_keys()) == set(EXPECTED)




# - - - BUG-6 / BUG-7 regression guard - - -
#
# These two assert the orientation fixed in Phase 4 on their own, so the
# change is covered by something that names it rather than only by a row in
# the table above. Both fail on every commit before the Phase 4 fix.

@pytest.mark.parametrize('key', ['DIFFERENTIAL SCANNING CALORIMETRY', 'LC/MS'])
def test_thermal_and_chromatographic_axes_run_forward(key, render, tmp_path):
    drawn = render(_mapping_keys()[key][0], tmp_path)
    low, high = drawn['xlim']
    assert low < high, (
        '%s must be drawn low -> high; reversed is the NMR convention and '
        'neither a thermal curve nor a chromatogram follows it' % key
    )


def test_dsc_matches_tga_orientation(render, tmp_path):
    """The defect was DSC diverging from its own family."""
    dsc = render(_mapping_keys()['DIFFERENTIAL SCANNING CALORIMETRY'][0], tmp_path)
    tga = render(_mapping_keys()['THERMOGRAVIMETRIC ANALYSIS'][0], tmp_path)
    assert dsc['orientation'] == tga['orientation'] == 'forward'


def test_unrecognised_datatype_rendering(render, tmp_path):
    """The generic curve path, which #291 established.

    A datatype absent from data_type.json is not an error: it renders with
    the default threshold, the reversed fallback orientation and generic
    axis labels.
    """
    drawn = render('NEUTRON SCATTERING', tmp_path)
    assert drawn['threshold'] == 0.5
    assert drawn['orientation'] == 'reversed'
    assert _x_style(drawn['xlabel']) == 'generic'
    assert _y_style(drawn['ylabel']) == 'generic'

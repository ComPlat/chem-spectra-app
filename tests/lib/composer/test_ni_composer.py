import pytest

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
# NIComposer selects the Agg backend on import, so import pyplot after it
from chem_spectra.lib.composer.ni import NIComposer
import matplotlib.pyplot as plt  # noqa: E402

source_nmr = './tests/fixtures/source/1H.dx'


@pytest.fixture
def ni_composer():
    base_converter = JcampBaseConverter(source_nmr)
    ni_converter = JcampNIConverter(base_converter)
    return NIComposer(ni_converter)


def test_ensure_itg_mpy_hydrates_from_core_tables(ni_composer):
    ni_composer.all_itgs = []
    ni_composer.mpys = []
    ni_composer.refShift = 0
    ni_composer.core.itg_table = ['(1.0, 2.0, 3.0)\n(4.0, 5.0, 6.0)\n']
    ni_composer.core.mpy_itg_table = ['(1, 1.0, 2.0, 0.0, 1.0, 0.0, m, 1.0)\n']
    ni_composer.core.mpy_pks_table = ['(1, 1.5, 100.0)\n(1, 1.7, 90.0)\n']

    ni_composer._ensure_itg_mpy_from_core_tables()

    assert ni_composer.all_itgs == [
        {'xL': 1.0, 'xU': 2.0, 'area': 3.0},
        {'xL': 4.0, 'xU': 5.0, 'area': 6.0},
    ]
    assert len(ni_composer.mpys) == 1
    mpy = ni_composer.mpys[0]
    assert mpy['mpyType'] == 'm'
    assert mpy['xExtent'] == {'xL': 1.0, 'xU': 2.0}
    assert mpy['peaks'] == [{'x': 1.5, 'y': 100.0}, {'x': 1.7, 'y': 90.0}]


def test_ensure_itg_mpy_skips_malformed_lines(ni_composer):
    ni_composer.all_itgs = []
    ni_composer.mpys = []
    ni_composer.core.itg_table = ['\n(1.0)\n(1.0, 2.0, 3.0)\n\n']
    ni_composer.core.mpy_itg_table = ['\n(1, 1.0, 2.0)\n(1, 1.0, 2.0, 0.0, 1.0, 0.0, s, 1.0)\n']
    ni_composer.core.mpy_pks_table = ['\n(1, 1.5, 100.0)\n(2)\n']

    ni_composer._ensure_itg_mpy_from_core_tables()

    assert ni_composer.all_itgs == [{'xL': 1.0, 'xU': 2.0, 'area': 3.0}]
    assert len(ni_composer.mpys) == 1
    assert ni_composer.mpys[0]['peaks'] == [{'x': 1.5, 'y': 100.0}]


def test_ensure_itg_mpy_skips_multiplet_without_peaks(ni_composer):
    # a multiplet whose index has no matching peaks must be dropped:
    # calc_mpy_center divides by the peak count when it is drawn later
    ni_composer.all_itgs = []
    ni_composer.mpys = []
    ni_composer.core.itg_table = []
    ni_composer.core.mpy_itg_table = ['(2, 1.0, 2.0, 0.0, 1.0, 0.0, t, 1.0)\n']
    ni_composer.core.mpy_pks_table = ['(1, 1.5, 100.0)\n']

    ni_composer._ensure_itg_mpy_from_core_tables()

    assert ni_composer.mpys == []


def test_ensure_itg_mpy_handles_missing_peaks_table(ni_composer):
    # $OBSERVEDMULTIPLETS present but $OBSERVEDMULTIPLETSPEAKS missing
    ni_composer.all_itgs = []
    ni_composer.mpys = []
    ni_composer.core.itg_table = []
    ni_composer.core.mpy_itg_table = ['(1, 1.0, 2.0, 0.0, 1.0, 0.0, m, 1.0)\n']
    ni_composer.core.mpy_pks_table = []

    ni_composer._ensure_itg_mpy_from_core_tables()

    assert ni_composer.mpys == []


def test_plot_overlays_returns_boundaries(ni_composer):
    ni_composer.core.itg_table = ['(1.0, 2.0, 3.0)\n']
    ni_composer.core.mpy_itg_table = ['(1, 1.0, 2.0, 0.0, 1.0, 0.0, m, 1.0)\n']
    ni_composer.core.mpy_pks_table = ['(1, 1.5, 100.0)\n']

    fig = plt.figure()
    try:
        y_boundary_min, y_boundary_max = ni_composer.plot_overlays(
            plt, adjust_xlim=False,
        )
    finally:
        plt.close(fig)

    assert y_boundary_min < y_boundary_max

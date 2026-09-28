"""The three polarity fields, each driven by data that actually exercises it.

`transmittance`, `peaks_inverted` and `negative_peaks` coincide for infrared
but express different concerns. Two of them failed nothing when first
introduced, because no fixture met their preconditions:

- `transmittance` inverts only when the median sits below half the maximum.
  `IR.dx` does not -- a real transmittance trace sits near 100%T with dips.
  So the probe here is an NMR-shaped file relabelled as infrared.
- `negative_peaks` only matters when the trace dips well below zero
  (`max_y * 0.4 < -min_y`). `13C-DEPT135.dx` is the only fixture in the repo
  that does.

Per the plan's rule a zero sensitivity count is a gap, not a pass.
"""

import re
from dataclasses import replace

import numpy as np

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

# median well below half the maximum: sharp peaks over a low baseline
ABSORBANCE_SHAPED = './tests/fixtures/source/1H.dx'
# the only fixture whose trace dips below zero
DIPPING = './tests/fixtures/source/13C-DEPT135.dx'


def _relabelled(source, datatype, tmp_path):
    body = open(source).read()
    header = re.search(r'##DATA TYPE=.*', body).group(0)
    target = tmp_path / 'probe.dx'
    target.write_text(body.replace(header, '##DATA TYPE=' + datatype, 1))
    return JcampTechniqueConverter(JcampBaseConverter(str(target)))


# - - - transmittance - - -

def test_transmittance_inverts_an_absorbance_shaped_trace(tmp_path):
    """Infrared and Raman differ only in this field, so ys differ only by it.

    Both are em-wave, so any orientation handling applies identically and
    cannot account for the difference.
    """
    infrared = _relabelled(ABSORBANCE_SHAPED, 'INFRARED SPECTRUM', tmp_path)
    raman = _relabelled(ABSORBANCE_SHAPED, 'RAMAN SPECTRUM', tmp_path)

    assert infrared.technique.transmittance is True
    assert raman.technique.transmittance is False

    raw = np.asarray(raman.ys, dtype=float)
    assert np.median(raw) < 0.5 * np.max(raw), 'probe no longer meets the precondition'
    assert np.allclose(np.asarray(infrared.ys, dtype=float), np.max(raw) - raw)


def test_transmittance_leaves_a_trace_that_is_already_transmittance(tmp_path):
    """The real IR fixture sits near 100%T, so nothing is inverted."""
    infrared = _relabelled('./tests/fixtures/source/IR.dx', 'INFRARED SPECTRUM', tmp_path)
    raw = np.asarray(infrared.ys, dtype=float)
    assert not np.median(raw) < 0.5 * np.max(raw)


# - - - negative_peaks - - -

def _peak_count(converter, negative_peaks):
    """Re-run peak picking with one field changed and nothing else."""
    converter.technique = replace(converter.technique, negative_peaks=negative_peaks)
    return len(converter._JcampTechniqueConverter__exec_peak_picking_logic())


def test_negative_peaks_folds_in_the_inverted_series(tmp_path):
    """A DEPT trace dips below zero; those peaks are real and must be kept."""
    converter = JcampTechniqueConverter(JcampBaseConverter(DIPPING))
    ys = np.asarray(converter.ys, dtype=float)
    assert np.max(ys) * 0.4 < -np.min(ys), 'probe no longer dips below zero'

    with_augmentation = _peak_count(converter, True)
    without = _peak_count(converter, False)
    assert with_augmentation > without, (
        'folding in the inverted series must find peaks the plain pass misses')


def test_infrared_and_cds_opt_out_of_the_augmentation():
    """Both opt out, for different reasons -- see the field comment.

    IR because peaks_inverted already searched the inverted series; CDS
    because its signal is genuinely bipolar and both lobes are already real
    peaks rather than artefacts to be recovered.
    """
    from chem_spectra.lib.converter.jcamp.techniques import SPECTRUM_TECHNIQUES
    assert SPECTRUM_TECHNIQUES['INFRARED'].negative_peaks is False
    assert SPECTRUM_TECHNIQUES['CIRCULAR DICHROISM SPECTROSCOPY'].negative_peaks is False
    assert SPECTRUM_TECHNIQUES['NMR'].negative_peaks is True


def test_cds_opt_out_is_load_bearing(tmp_path):
    """Circular dichroism really would gain peaks without the opt-out."""
    converter = _relabelled(DIPPING, 'CIRCULAR DICHROISM SPECTROSCOPY', tmp_path)
    assert converter.technique.negative_peaks is False
    baseline = len(converter._JcampTechniqueConverter__exec_peak_picking_logic())
    assert _peak_count(converter, True) > baseline


def test_infrared_opt_out_is_redundant_not_load_bearing(tmp_path):
    """Measured, and worth stating: for infrared the opt-out changes nothing.

    `peaks_inverted` already runs find_peaks on 1 - ys, so the augmentation
    has nothing left to discover and flipping `negative_peaks` back on gives
    an identical peak set. The field is set False for infrared to record the
    intent, not because the behaviour depends on it.

    If this test ever starts failing, the two fields have stopped
    overlapping -- which is the case the field comment warns about, and means
    infrared's value must be decided on its own merits rather than inherited.
    """
    converter = _relabelled(DIPPING, 'INFRARED SPECTRUM', tmp_path)
    assert converter.technique.negative_peaks is False
    assert converter.technique.peaks_inverted is True
    baseline = len(converter._JcampTechniqueConverter__exec_peak_picking_logic())
    assert _peak_count(converter, True) == baseline

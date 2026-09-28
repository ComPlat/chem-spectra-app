"""A sorption isotherm's two branches are labelled and marked by direction.

`tf_combine` (and the bagit writer, identically) renames a sorption trace to
ADSORPTION or DESORPTION and marks it '^' or 'v' depending on which way x
runs. No fixture declares the technique, so the branch was never executed by
the suite and flipping `sorption_branches` off failed nothing -- a gap, not
a pass.

The probe is `IR.dx` relabelled -- it stores x descending, and exchanging
its declared ends makes it ascending. It must not be an `(XY..XY)` /
XYPOINTS file: for those `__read_xs` takes x from the data table and ignores
FIRSTX/LASTX entirely, so the direction cannot be varied that way.
"""

import io
import re

from werkzeug.datastructures import FileStorage

import chem_spectra.model.transformer as transformer_module
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel

SOURCE = './tests/fixtures/source/IR.dx'   # XYDATA: honours FIRSTX/LASTX
MOLFILE = './tests/fixtures/source/molfile/svs813f1_B.mol'


def _probe_bytes(datatype, descending=True):
    """IR.dx stores x descending; ascending is the exchanged variant."""
    body = open(SOURCE).read()
    body = body.replace(
        re.search(r'##DATA TYPE=.*', body).group(0), '##DATA TYPE=' + datatype, 1)
    if not descending:
        body = body.replace('##FIRSTX=3997.453', '##FIRSTX=@@') \
                   .replace('##LASTX=373.96442', '##LASTX=3997.453') \
                   .replace('##FIRSTX=@@', '##FIRSTX=373.96442') \
                   .replace('##DELTAX=-1.4165319', '##DELTAX=1.4165319')
    return body.encode()


def _plotted(datatype, descending, monkeypatch):
    """Run tf_combine and report the label/marker each trace was drawn with."""
    drawn = []
    real_plot = transformer_module.plt.plot

    def spy_plot(*args, **kwargs):
        if 'label' in kwargs:
            drawn.append((kwargs.get('label'), kwargs.get('marker')))
        return real_plot(*args, **kwargs)

    monkeypatch.setattr(transformer_module.plt, 'plot', spy_plot)

    files = [
        FileContainer(FileStorage(io.BytesIO(_probe_bytes(datatype, descending)),
                                  filename='probe.dx'))
        for _ in range(2)
    ]
    with open(MOLFILE) as molfile:
        TransformerModel(None, molfile=molfile, params={'ext': 'jdx'},
                         multiple_files=files).tf_combine()
    return drawn


def test_ascending_sorption_trace_is_labelled_adsorption(monkeypatch):
    drawn = _plotted('SORPTION-DESORPTION MEASUREMENT', False, monkeypatch)
    assert drawn, 'nothing was plotted'
    assert all(label == 'ADSORPTION' and marker == '^' for label, marker in drawn)


def test_descending_sorption_trace_is_labelled_desorption(monkeypatch):
    drawn = _plotted('SORPTION-DESORPTION MEASUREMENT', True, monkeypatch)
    assert drawn, 'nothing was plotted'
    assert all(label == 'DESORPTION' and marker == 'v' for label, marker in drawn)


def test_another_technique_keeps_its_filename_and_no_marker(monkeypatch):
    """The control: same bytes, same direction, a technique without branches."""
    drawn = _plotted('HPLC UV/VIS SPECTRUM', False, monkeypatch)
    assert drawn, 'nothing was plotted'
    assert all(label not in ('ADSORPTION', 'DESORPTION') for label, _ in drawn)
    assert all(marker == '' for _, marker in drawn)

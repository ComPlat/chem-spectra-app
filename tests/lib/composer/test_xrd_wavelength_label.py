"""The XRD x-axis label carries the wavelength's own unit.

It was `"WL={} nm".format(value, unit)`: the unit was passed and ignored, so
a wavelength given in Å was labelled in nm.
"""

import json
import shutil

import pytest

import chem_spectra.lib.composer.technique as technique_module
import chem_spectra.lib.converter.bagit.base as bagit_module
from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.bagit.base import BagItBaseConverter
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

XRD = './tests/fixtures/source/xrd/Test_data_op_001_XRD.jdx'
ANGSTROM = {'name': 'CuKalpha', 'value': 1.5406, 'label': 'Cu K-alpha', 'unit': 'Å'}


def _xlabels(module, monkeypatch):
    labels = []
    real = module.plt.xlabel

    def spy(text, **kwargs):
        labels.append(text)
        return real(text, **kwargs)

    monkeypatch.setattr(module.plt, 'xlabel', spy)
    return labels


@pytest.mark.parametrize('wavelength, expected', [
    (None, 'WL=0.15406 nm'),
    (ANGSTROM, 'WL=1.5406 Å'),
])
def test_the_single_image_label(monkeypatch, wavelength, expected):
    params = {'waveLength': json.dumps(wavelength)} if wavelength else {}
    labels = _xlabels(technique_module, monkeypatch)
    composer = TechniqueComposer(
        JcampTechniqueConverter(JcampBaseConverter(XRD, params)))
    composer.tf_img().close()
    assert labels and labels[-1].endswith(expected)


def test_the_archive_overlay_label(tmp_path, monkeypatch):
    data = tmp_path / 'data'
    data.mkdir()
    for name in ('a.jdx', 'b.jdx'):
        shutil.copy(XRD, data / name)
    labels = _xlabels(bagit_module, monkeypatch)
    BagItBaseConverter(str(tmp_path), {'waveLength': json.dumps(ANGSTROM)})
    assert labels and labels[-1].endswith('WL=1.5406 Å')

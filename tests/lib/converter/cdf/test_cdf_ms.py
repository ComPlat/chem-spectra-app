"""An ANDI/MS netCDF upload keeps every scan, the last one included.

`scan_index` holds each scan's first point. The runs were built from pairs
of consecutive scan starts, so the last scan, which has no successor, was
dropped from every ANDI/MS upload. The netCDF here is synthetic.
"""
import netCDF4
import numpy as np
from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.model.transformer import TransformerModel


def _andi_ms_cdf(path, scans):
    """A minimal ANDI/MS netCDF: `scans` is a list of (mzs, intensities)."""
    with netCDF4.Dataset(path, 'w') as cdf:
        cdf.dataset_completeness = 'C1+C2'
        cdf.ms_template_revision = '1.0.1'
        n_points = sum(len(mzs) for mzs, _ in scans)
        cdf.createDimension('point_number', n_points)
        cdf.createDimension('scan_number', len(scans))
        mass = cdf.createVariable('mass_values', 'f8', ('point_number',))
        intensity = cdf.createVariable(
            'intensity_values', 'f4', ('point_number',))
        index = cdf.createVariable('scan_index', 'i4', ('scan_number',))
        mass[:] = np.concatenate([mzs for mzs, _ in scans])
        intensity[:] = np.concatenate([ys for _, ys in scans])
        index[:] = np.cumsum([0] + [len(mzs) for mzs, _ in scans[:-1]])


def test_every_scan_arrives_the_last_one_included(tmp_path):
    scans = [
        (np.array([100.0, 101.0, 102.0]), np.array([5.0, 50.0, 5.0])),
        (np.array([200.0, 201.0]), np.array([7.0, 70.0])),
        (np.array([300.0, 301.0, 302.0, 303.0]),
         np.array([1.0, 9.0, 90.0, 1.0])),
    ]
    path = tmp_path / 'ms.cdf'
    _andi_ms_cdf(path, scans)

    with open(path, 'rb') as f:
        file = FileContainer(FileStorage(f, filename='ms.cdf'))
        composer, invalid_molfile = TransformerModel(
            file, params={'fname': 'ms.cdf', 'ext': 'cdf'}).to_composer()

    assert isinstance(composer, MSComposer)
    assert invalid_molfile is False
    spectra = composer.core.spectra
    assert [len(s) for s in spectra] == [3, 2, 4]
    assert spectra[-1][:, 0].tolist() == [300.0, 301.0, 302.0, 303.0]
    jcamp = composer.tf_jcamp().read().decode('utf-8').split('\n')
    assert '##$CSSCANCOUNT=3' in jcamp

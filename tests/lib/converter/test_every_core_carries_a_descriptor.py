"""Every converter a composer is given carries a spectrum-technique descriptor.

The composer decides integrals, multiplicity, NMR headers and axis direction
from `core.technique`. A core without one is drawn as the generic curve, and
nothing raises: the NMRium converter went without one for years, resolved to
NMR only by a fallback on its `non_nmr` flag.

Every `*Converter` class under `chem_spectra/lib/converter` is listed below,
either as a composer core (constructed here, its descriptor checked) or as
not one, with the reason. A converter added later fails
`test_every_converter_class_is_accounted_for` until it is placed in one list.
"""
import importlib
import inspect
import pkgutil

import netCDF4
import numpy as np
import pytest
from werkzeug.datastructures import FileStorage

import chem_spectra.lib.converter as converter_package
from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.lib.converter.cdf.base import CdfBaseConverter
from chem_spectra.lib.converter.cdf.ms import CdfMSConverter
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.jcamp.curve import JcampCurveConverter
from chem_spectra.lib.converter.jcamp.techniques import (
    SPECTRUM_TECHNIQUES, SpectrumTechnique,
)
from chem_spectra.lib.converter.ms import MSConverter
from chem_spectra.lib.converter.nmrium.base import NMRiumDataConverter

FIXTURES = './tests/fixtures/source/'
MS_PARAMS = {'mass': 230.079907196}


def _jcamp_technique(_tmp_path):
    return JcampCurveConverter(JcampBaseConverter(FIXTURES + '1H.dx'))


def _jcamp_ms(_tmp_path):
    return JcampMSConverter(
        JcampBaseConverter(FIXTURES + 'ms/svs813f1.jdx', MS_PARAMS))


def _mzml(_tmp_path):
    with open(FIXTURES + 'ms/svs813f1.mzML', 'rb') as f:
        return MSConverter(FileContainer(FileStorage(f)), MS_PARAMS)


def _andi_ms(tmp_path):
    path = tmp_path / 'ms.cdf'
    with netCDF4.Dataset(path, 'w') as cdf:
        cdf.createDimension('point_number', 4)
        cdf.createDimension('scan_number', 2)
        mass = cdf.createVariable('mass_values', 'f8', ('point_number',))
        intensity = cdf.createVariable(
            'intensity_values', 'f4', ('point_number',))
        index = cdf.createVariable('scan_index', 'i4', ('scan_number',))
        mass[:] = np.array([100.0, 101.0, 200.0, 201.0])
        intensity[:] = np.array([5.0, 50.0, 7.0, 70.0])
        index[:] = np.array([0, 2])
    return CdfMSConverter(CdfBaseConverter(str(path), {'fname': 'ms.cdf'}))


def _nmrium(_tmp_path):
    with open(FIXTURES + 'nmrium/flat.nmrium', 'rb') as f:
        return NMRiumDataConverter(
            FileContainer(FileStorage(f, filename='spectrum.nmrium')))


# The converters handed to CurveComposer or MSComposer, with what each
# must resolve to for the fixture it is built from.
CORES = {
    'JcampCurveConverter': (_jcamp_technique, 'NMR'),
    'JcampMSConverter': (_jcamp_ms, 'MS'),
    'MSConverter': (_mzml, 'MS'),
    'CdfMSConverter': (_andi_ms, 'MS'),
    'NMRiumDataConverter': (_nmrium, 'NMR'),
}

# Converters that never reach a composer themselves.
NOT_CORES = {
    'JcampBaseConverter': 'wrapped by JcampCurveConverter or JcampMSConverter',
    'FidBaseConverter': 'wrapped by JcampCurveConverter, which reads its descriptor',
    'CdfBaseConverter': 'wrapped by CdfMSConverter',
    'BagItBaseConverter': 'an archive of cores, not a core',
}


def _converter_classes():
    found = set()
    for module_info in pkgutil.walk_packages(
            converter_package.__path__, converter_package.__name__ + '.'):
        module = importlib.import_module(module_info.name)
        for name, obj in inspect.getmembers(module, inspect.isclass):
            if obj.__module__ == module.__name__ and name.endswith('Converter'):
                found.add(name)
    return found


def test_every_converter_class_is_accounted_for():
    assert _converter_classes() == set(CORES) | set(NOT_CORES)


@pytest.mark.parametrize('name', sorted(CORES))
def test_each_core_carries_its_descriptor(name, tmp_path):
    build, key = CORES[name]

    core = build(tmp_path)

    assert isinstance(core.technique, SpectrumTechnique)
    assert core.technique is SPECTRUM_TECHNIQUES[key]

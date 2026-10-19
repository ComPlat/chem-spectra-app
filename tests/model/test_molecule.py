"""The molfile's temp file is closed once RDKit has read it.

`tf.close` without the call left it open until garbage collection.
"""

import chem_spectra.model.molecule as molecule_module
from chem_spectra.model.molecule import MoleculeModel

MOLFILE = './tests/fixtures/source/molfile/svs813f1_B.mol'


def test_the_temp_molfile_is_closed(monkeypatch):
    opened = []
    real = molecule_module.store_str_in_tmp

    def spy(*args, **kwargs):
        tf = real(*args, **kwargs)
        opened.append(tf)
        return tf

    monkeypatch.setattr(molecule_module, 'store_str_in_tmp', spy)
    with open(MOLFILE) as handle:
        model = MoleculeModel(handle.read())

    assert model.mol is not None
    assert opened and all(tf.closed for tf in opened)

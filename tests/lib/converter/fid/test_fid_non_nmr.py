"""`FidBaseConverter` states NMR literally rather than deriving it.

A free induction decay is an NMR measurement by definition, and the
converter hardcodes `typ = 'NMR'`. It used to derive `non_nmr` from a
predicate testing that hardcoded value against ten *other* techniques,
which could never match -- dead code that read like live dispatch, and
stale by seven data_type.json keys besides.
"""

from chem_spectra.lib.converter.fid.base import FidBaseConverter


def test_fid_converter_states_nmr_literally():
    source = open('chem_spectra/lib/converter/fid/base.py', encoding='utf-8').read()
    assert "self.typ = 'NMR'" in source
    assert "self.datatypes = ['NMR SPECTRUM']" in source
    assert 'self.non_nmr = False' in source


def test_fid_has_no_unreachable_non_nmr_predicate():
    """The deleted predicate must not come back.

    Anything re-deriving non_nmr from a technique list here is either dead
    (typ is hardcoded) or a sign that typ became dynamic without the list
    being completed.
    """
    assert not hasattr(FidBaseConverter, '_FidBaseConverter__non_nmr')

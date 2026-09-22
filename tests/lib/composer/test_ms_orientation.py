"""A mass spectrum renders ascending, so the registry must say so.

`MS`'s `x_reversed` is already `False` -- corrected in #293 (`3bbdc2f`)
after the rendering matrix had measured it from a path MS never takes. It is
still unread, because nothing routes MS through `TechniqueComposer`, the
field's only consumer.

So this pins a value rather than fixing one, and it pins it against the
measurement instead of against the table: `MSComposer.tf_img` never calls
`plt.xlim`, so matplotlib autoscales, and all three fixtures render
ascending -- m/z low to high, as a mass spectrum is conventionally drawn.
Routing MS through the shared composer makes the field live, and a pin
measured from a render is what stops that silently mirroring every mass
spectrum.

Note the axis limits have to be read *during* `savefig`: `tf_img` calls
`plt.clf()` immediately after, so reading them afterwards returns the
default (0, 1) and proves nothing.
"""

import pytest

import chem_spectra.lib.composer.ms as ms_module
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.converter.jcamp.techniques import SPECTRUM_TECHNIQUES

MS_FIXTURES = [
    './tests/fixtures/source/ms/svs813f1.jdx',
    './tests/fixtures/source/ms/MS_ESI.jdx',
    './tests/fixtures/source/ms/ms_v6.dx',
]


def _rendered_xlim(composer):
    captured = {}
    real_savefig = ms_module.plt.savefig

    def spy_savefig(*args, **kwargs):
        captured['xlim'] = ms_module.plt.gca().get_xlim()
        return real_savefig(*args, **kwargs)

    ms_module.plt.savefig = spy_savefig
    try:
        composer.tf_img().close()
    finally:
        ms_module.plt.savefig = real_savefig
    return captured['xlim']


def test_ms_registry_entry_is_not_reversed():
    assert SPECTRUM_TECHNIQUES['MS'].x_reversed is False


@pytest.mark.parametrize('fixture', MS_FIXTURES)
def test_an_ms_render_is_ascending(fixture):
    """The measurement the registry value has to agree with."""
    composer = MSComposer(JcampMSConverter(JcampBaseConverter(fixture)))
    low, high = _rendered_xlim(composer)
    assert low < high, '%s rendered descending' % fixture


def test_ms_composer_never_sets_xlim():
    """Why the render is ascending: nothing constrains the axis.

    If this ever starts failing, MSComposer has begun choosing an
    orientation of its own and the registry value has to be re-derived from
    whatever it chooses.
    """
    called = []
    real_xlim = ms_module.plt.xlim

    def spy_xlim(*args, **kwargs):
        if args:
            called.append(args)
        return real_xlim(*args, **kwargs)

    ms_module.plt.xlim = spy_xlim
    try:
        composer = MSComposer(JcampMSConverter(JcampBaseConverter(MS_FIXTURES[0])))
        composer.tf_img().close()
    finally:
        ms_module.plt.xlim = real_xlim
    assert called == []

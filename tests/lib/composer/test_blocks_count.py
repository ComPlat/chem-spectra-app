"""`##BLOCKS` states how many child blocks the LINK actually contains.

It was hardcoded to 1 with a `# TBD`, while a 1H file carries three, so every
file this app produced misdeclared its own structure.

The count cannot be known when the root header is written -- the children are
generated afterwards -- so the body is composed first and the header prepended
once the total is known. It is derived by counting `##END=` lines rather than
tracked as the body is built, so a branch added later cannot forget it.

**Not a fix for issue #290.** That reports a blank chart and blames this
record; neither resolved `jcampconverter` version reads `##BLOCKS` at all, and
PR #289 already retracted that claim once. This is correctness on its own.
"""

import re

import pytest

from chem_spectra.lib.composer.base import BaseComposer
from chem_spectra.lib.composer.technique import TechniqueComposer
from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter

SOURCES = [
    './tests/fixtures/source/1H.dx',
    './tests/fixtures/source/13C-DEPT135.dx',
    './tests/fixtures/source/IR.dx',
    './tests/fixtures/source/JPK-948.jdx',
    './tests/fixtures/source/hplc/chromatogram.jdx',
]


def _meta(source):
    return ''.join(TechniqueComposer(
        JcampTechniqueConverter(JcampBaseConverter(source))).meta)


@pytest.mark.parametrize('source', SOURCES)
def test_declared_blocks_matches_the_children_written(source):
    meta = _meta(source)
    declared = int(re.search(r'^##BLOCKS=(\d+)', meta, re.M).group(1))
    endings = len(re.findall(r'^##END=', meta, re.M))
    # one ##END= closes the LINK itself; the rest close children
    assert declared == endings - 1


@pytest.mark.parametrize('source', SOURCES)
def test_it_is_no_longer_always_one(source):
    """The specific defect: every layout declared 1 regardless of content."""
    meta = _meta(source)
    assert int(re.search(r'^##BLOCKS=(\d+)', meta, re.M).group(1)) > 1


def test_the_count_is_derived_from_the_body():
    """Deriving it is what makes it survive a branch added later."""
    body = ['##TITLE=x\n', '##END=\n', '##TITLE=y\n', '##END=\n', '##END=\n']
    assert BaseComposer.count_child_blocks(body) == 2


def test_a_body_with_no_children_still_declares_one():
    """Degenerate input must not emit ##BLOCKS=0 or a negative number."""
    assert BaseComposer.count_child_blocks(['##END=\n']) == 1
    assert BaseComposer.count_child_blocks([]) == 1

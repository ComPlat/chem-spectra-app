"""Legend labels in the combined image are bounded (issue #265).

Long filenames ran across the plot and hid the curves. The middle is dropped
rather than the tail, because names in this domain often differ only by their
suffix -- `..._13C.jdx` against `..._1H.jdx` -- so truncating the end alone
can make two curves indistinguishable in the legend, which is the problem the
issue is about.
"""

import pytest

from chem_spectra.lib.shared.misc import LEGEND_LABEL_MAX, shorten_label


@pytest.mark.parametrize('label', [
    'short.jdx',
    'SVS-790A_13C.jdx',
    'x' * LEGEND_LABEL_MAX,
])
def test_labels_within_the_limit_are_untouched(label):
    assert shorten_label(label) == label


def test_a_long_label_is_bounded():
    assert len(shorten_label('y' * 200)) == LEGEND_LABEL_MAX


def test_both_ends_survive():
    """The discriminating part of these names is the suffix."""
    a = shorten_label('a-very-long-sample-identifier-from-the-instrument_13C.jdx')
    b = shorten_label('a-very-long-sample-identifier-from-the-instrument_1H.jdx')
    assert a != b, 'two curves would be indistinguishable in the legend'
    assert a.endswith('.jdx') and b.endswith('.jdx')
    assert a.startswith('a-very')


def test_every_fixture_name_is_unaffected():
    """The limit is chosen so the common case never changes."""
    import glob
    import os
    names = [os.path.basename(p)
             for p in glob.glob('tests/fixtures/**/*.*', recursive=True)]
    assert names
    assert all(shorten_label(n) == n for n in names)


def test_none_passes_through():
    assert shorten_label(None) is None


def test_the_combined_image_actually_uses_it(tmp_path):
    """Asserted through tf_combine, not just on the helper.

    The bug is about what reaches `plt.plot(label=...)`; a helper-only test
    would pass with the call site unchanged.
    """
    import io
    import chem_spectra.model.transformer as transformer_module
    from werkzeug.datastructures import FileStorage
    from chem_spectra.controller.helper.file_container import FileContainer
    from chem_spectra.model.transformer import TransformerModel

    drawn = []
    real_plot = transformer_module.plt.plot

    def spy_plot(*args, **kwargs):
        if 'label' in kwargs:
            drawn.append(kwargs['label'])
        return real_plot(*args, **kwargs)

    long_name = 'a-very-long-sample-identifier-from-the-instrument-run-seven_13C.jdx'
    payload = open('./tests/fixtures/source/1H.dx', 'rb').read()
    files = [
        FileContainer(FileStorage(io.BytesIO(payload), filename=long_name))
        for _ in range(2)
    ]

    transformer_module.plt.plot = spy_plot
    try:
        with open('./tests/fixtures/source/molfile/svs813f1_B.mol') as molfile:
            TransformerModel(None, molfile=molfile, params={'ext': 'jdx'},
                             multiple_files=files).tf_combine()
    finally:
        transformer_module.plt.plot = real_plot

    assert drawn, 'nothing was plotted'
    assert all(len(label) <= LEGEND_LABEL_MAX for label in drawn)
    assert drawn[0].endswith('_13C.jdx')


def test_the_bagit_combined_image_uses_names_not_indices():
    """The BagIt legend read `0, 1, 2 ...` -- it never showed a filename.

    `__combine_images` took a `list_file_names` parameter that its only caller
    never passed, so the branch was dead. The name is now carried on the
    composer instead of in an index-aligned list, because that method filters
    LC/MS and MS composers out again and any parallel list would mislabel the
    survivors.
    """
    import tempfile
    import zipfile
    import chem_spectra.lib.converter.bagit.base as bagit_module
    from chem_spectra.lib.converter.bagit.base import BagItBaseConverter

    drawn = []
    real_plot = bagit_module.plt.plot

    def spy_plot(*args, **kwargs):
        if 'label' in kwargs:
            drawn.append(kwargs['label'])
        return real_plot(*args, **kwargs)

    bagit_module.plt.plot = spy_plot
    try:
        with tempfile.TemporaryDirectory() as td:
            archive = './tests/fixtures/source/bagit/cv/File053_BagIt.zip'
            with zipfile.ZipFile(archive, 'r') as z:
                z.extractall(td)
            BagItBaseConverter(td).combined_image
    finally:
        bagit_module.plt.plot = real_plot

    assert drawn, 'nothing was plotted'
    assert not all(label.isdigit() for label in drawn), (
        'legend fell back to indices; the filename is not reaching the plot')
    assert all(len(label) <= LEGEND_LABEL_MAX for label in drawn)

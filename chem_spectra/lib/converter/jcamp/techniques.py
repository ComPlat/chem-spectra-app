"""One descriptor per measurement technique.

`typ` is a single-valued fact -- which technique a file represents -- that
the codebase re-derives as sixteen independent booleans
(`is_xrd`, `is_cyclic_volta`, ... plus the `is_em_wave` grouping and the
`non_nmr` negation) and then re-assembles into if-chains at 128 call sites
across seven files. Adding a technique means editing five of them -- which
is how DSC and LC/MS ended up drawn mirrored (fixed in #292) and how
`data_type.json` drifted twelve entries behind chemotion-converter-app
(fixed in #291). Both were symptoms of this design.

This module holds what those chains decide, keyed by the `data_type.json`
key, so a technique becomes one mapping line plus one entry here.

The axis and threshold values were extracted from the pre-refactor tree by
execution, not by reading, and are pinned by
`tests/lib/composer/test_technique_rendering_matrix.py`.

The four NMR gates are consumed by `BaseComposer` and `TechniqueComposer`,
and their behaviour is covered by `tests/lib/composer/test_non_nmr_gates.py`
-- flipping any of them to the wrong value fails the suite. Two fields are
**not** yet read and so are pinned only by table-to-table parity in
`tests/lib/converter/jcamp/test_techniques.py`:

- `cv_scaling`, which nothing consumes at all;
  (`x_axis='xrd'` was in this state too until the composer's hardcoded
  `is_xrd` branch was migrated to read it -- see review of #293);
- `threshold`, because `converter/jcamp/technique.py` still carries its own
  threshold table.

Both need render-level assertions when they are wired up.

Note on the NMR gates: `non_nmr` gated four unrelated concerns --
integration pairing, multiplicity output, axis-label style and peak
annotation. They are separate fields here on purpose. Today all four are
True only for NMR, so collapsing them into one boolean would reproduce
current behaviour exactly while quietly asserting that they must always
agree, which is not something the code establishes.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class SpectrumTechnique:
    """How one technique is parsed and drawn."""

    key: str

    # - - - axis presentation - - -
    x_axis: str = 'generic'        # 'chemical_shift' | 'generic' | 'raw' | 'xrd'
    y_axis: str = 'generic'        # 'intensity'      | 'generic' | 'raw'
    # high -> low, the NMR convention. True is the fallback an unrecognised
    # datatype gets today, so it is the default here too.
    x_reversed: bool = True

    # - - - peak detection - - -
    threshold: float = 0.5

    # - - - the four concerns `non_nmr` gates, deliberately separate - - -
    # composer/base.py prepare_itg_mpy: pair integrations with multiplets
    nmr_integration: bool = False
    # composer/base.py gen_mpy_*_info, composer/technique.py meta and drawing
    multiplicity: bool = False
    # composer/technique.py __draw_peaks: annotate peaks on the plot
    peak_annotation: bool = False
    # composer/technique.py __gen_headers_spectrum_orig: the fuller NMR header
    nmr_headers: bool = False

    # - - - groupings and per-technique extras - - -
    # IR / Raman / UV-Vis share a header shape (the old `is_em_wave`)
    em_wave: bool = False
    # cyclic voltammetry y-scaling and the shared 10^n axis label
    cv_scaling: bool = False


def _nmr(key):
    return SpectrumTechnique(
        key,
        x_axis='chemical_shift', y_axis='intensity',
        x_reversed=True, threshold=0.005,
        nmr_integration=True, multiplicity=True,
        peak_annotation=True, nmr_headers=True,
    )


SPECTRUM_TECHNIQUES = {
    'NMR': _nmr('NMR'),

    'INFRARED': SpectrumTechnique('INFRARED', x_reversed=True, threshold=0.93,
                             em_wave=True),
    'RAMAN': SpectrumTechnique('RAMAN', x_reversed=True, threshold=0.07,
                          em_wave=True),
    # MS is not routed through TechniqueComposer yet: every production
    # `typ == 'MS'` path goes to MSComposer (transformer.py:273, :381,
    # bagit/base.py:82), which draws sticks and never calls plt.xlim, so
    # m/z renders ascending. This entry is therefore unread today, and
    # x_reversed=False is what it must be when the MS fold makes it live.
    'MS': SpectrumTechnique('MS', x_reversed=False, threshold=0.05),

    'HPLC UVVIS': SpectrumTechnique('HPLC UVVIS', x_reversed=False, threshold=0.05),
    # Also covers Avantes AvaSoft exports (data_type.json), whose instruments
    # are UV/VIS/NIR -- wider than this name, but handled identically: nm on x
    # running forward, same threshold. No NIR technique exists to route to.
    'UVVIS': SpectrumTechnique('UVVIS', x_reversed=False, threshold=0.05,
                          em_wave=True),

    'THERMOGRAVIMETRIC ANALYSIS': SpectrumTechnique(
        'THERMOGRAVIMETRIC ANALYSIS', x_reversed=False, threshold=1.05),
    'X-RAY DIFFRACTION': SpectrumTechnique(
        'X-RAY DIFFRACTION', x_axis='xrd', x_reversed=False, threshold=1.00),
    'CYCLIC VOLTAMMETRY': SpectrumTechnique(
        'CYCLIC VOLTAMMETRY', x_axis='raw', y_axis='raw',
        x_reversed=False, threshold=1.00, cv_scaling=True),
    'SIZE EXCLUSION CHROMATOGRAPHY': SpectrumTechnique(
        'SIZE EXCLUSION CHROMATOGRAPHY', x_reversed=False, threshold=0.5),
    'CIRCULAR DICHROISM SPECTROSCOPY': SpectrumTechnique(
        'CIRCULAR DICHROISM SPECTROSCOPY', x_reversed=False, threshold=1.00),
    'SORPTION-DESORPTION MEASUREMENT': SpectrumTechnique(
        'SORPTION-DESORPTION MEASUREMENT', x_reversed=False, threshold=1.00),
    'Emissions': SpectrumTechnique('Emissions', x_reversed=False, threshold=0.5),
    'DLS ACF': SpectrumTechnique('DLS ACF', x_reversed=False, threshold=1.05),
    'DLS intensity': SpectrumTechnique('DLS intensity', x_reversed=False,
                                  threshold=1.00),

    # Forward, like TGA. It was reversed until #292, because `is_dsc` was
    # never added to the hand-maintained orientation chain in tf_img.
    'DIFFERENTIAL SCANNING CALORIMETRY': SpectrumTechnique(
        'DIFFERENTIAL SCANNING CALORIMETRY', x_reversed=False, threshold=1.05),

    'GAS CHROMATOGRAPHY': SpectrumTechnique('GAS CHROMATOGRAPHY', x_reversed=False,
                                       threshold=0.5),

    # Forward: a chromatogram runs forward in retention time. It was
    # reversed until #292, because LC/MS has no is_* flag of its own and so
    # inherited the reversed default.
    'LC/MS': SpectrumTechnique('LC/MS', x_reversed=False, threshold=0.5),
}


# What an unrecognised datatype gets: the generic curve path established by
# b87ea91. Its key is '' so that `typ` and `technique.key` stay in step.
UNKNOWN_TECHNIQUE = SpectrumTechnique('')


def technique_for(typ):
    """The descriptor for a `typ`, falling back to the generic curve."""
    return SPECTRUM_TECHNIQUES.get(typ, UNKNOWN_TECHNIQUE)

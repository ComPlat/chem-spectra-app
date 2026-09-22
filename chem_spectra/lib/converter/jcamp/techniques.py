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
-- flipping any of them to the wrong value fails the suite.

`threshold` is read by `converter/jcamp/technique.py`, which used to carry a
duplicate table keyed on the raw datatype string. Flipping NMR's value fails
13 tests outside the parity file, so it is load-bearing on the peak-detection
path, not just pinned by parity.

`x_axis` was in the parity-only state too for `'xrd'`, until the composer's
hardcoded `is_xrd` branch was migrated to read it -- see the review of #293.

`cyclic_voltammetry` (formerly `cv_scaling`, which nothing consumed) is now
read by all fourteen CV sites. Every field in this dataclass is read by
something.

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
    # this technique is cyclic voltammetry. Fourteen sites across the
    # converter, both composers, the bagit writer and the transformer ask
    # exactly that question -- the shift offset, the display info, the data
    # table, the y-scaling and the axis label. It was named `cv_scaling` for
    # only the last of those and went unconsumed; the name now matches what
    # it is asked.
    cyclic_voltammetry: bool = False

    # - - - signal polarity, three concerns that coincide for infrared - - -
    # converter/jcamp/technique.py __read_ys: a transmittance trace stored the
    # absorbance way up is inverted, judged by median against max.
    transmittance: bool = False
    # converter/jcamp/technique.py __exec_peak_picking_logic and
    # __run_auto_pick_peak: bands are troughs, so find_peaks runs on 1 - ys
    # and the auto table keeps the *lowest* hundred.
    peaks_inverted: bool = False
    # converter/jcamp/technique.py __exec_peak_picking_logic: fold in peaks
    # found on the inverted series when the trace dips well below zero -- the
    # DEPT case. Infrared opts out because peaks_inverted already did that
    # work; circular dichroism opts out because its signal is genuinely
    # bipolar and both lobes are real. Same effect, different reasons, so if
    # one of them ever changes this field is the wrong place to express it.
    negative_peaks: bool = True
    # converter/jcamp/technique.py __set_label: a y-axis declaring absorbance
    # is reported as absorbance rather than rewritten to TRANSMITTANCE.
    # True only for UVVIS, which reproduces the pre-refactor `not is_uv_vis`
    # guard exactly. NOTE: 'HPLC UVVIS' is a separate key and so does *not*
    # get this, meaning an HPLC file declaring ##YUNITS=ABSORBANCE is
    # relabelled TRANSMITTANCE -- the inverse quantity. Reachable, pinned by
    # test_absorbance_label.py, and left as-is: it is a domain call, not the
    # refactor's to make. See CHANGELOG.refactor-finish-flag-migration.md.
    absorbance_label: bool = False

    # - - - the two UV/VIS concerns, which cover different sets - - -
    # composer/base.py _build_integration_lines and prepare_itg_mpy: the
    # integration table gains an AUC column. HPLC UV/VIS only.
    auc_column: bool = False
    # composer/base.py _supports_visual_split and composer/technique.py
    # __uses_auc_drawing: visual integration splits, and integrations drawn
    # as areas. HPLC UV/VIS *and* plain UV/VIS -- a wider set than
    # auc_column, which is why these are two fields.
    visual_split: bool = False

    # - - - size exclusion chromatography and differential scanning - - -
    # composer/technique.py __generate_info_box: which annotation box is
    # drawn on the plot. '' draws none; 'sec' lists MN/MW/MP/D, 'dsc' lists
    # the melting point and Tg. A discriminator rather than two booleans,
    # because the box has exactly one variant per technique.
    info_box: str = ''
    # composer/technique.py __gen_header_sec: the SEC header block
    sec_headers: bool = False
    # composer/technique.py __gen_header_user_input_meta_data: melting point
    # and Tg written into the JCAMP, from params or from the file's own LDRs
    dsc_metadata: bool = False


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
                             em_wave=True, transmittance=True,
                             peaks_inverted=True, negative_peaks=False),
    'RAMAN': SpectrumTechnique('RAMAN', x_reversed=True, threshold=0.07,
                          em_wave=True),
    # MS is not routed through TechniqueComposer yet: every production
    # `typ == 'MS'` path goes to MSComposer (transformer.py:273, :381,
    # bagit/base.py:82), which draws sticks and never calls plt.xlim, so
    # m/z renders ascending. This entry is therefore unread today, and
    # x_reversed=False is what it must be when the MS fold makes it live.
    'MS': SpectrumTechnique('MS', x_reversed=False, threshold=0.05),

    'HPLC UVVIS': SpectrumTechnique('HPLC UVVIS', x_reversed=False, threshold=0.05,
                               auc_column=True, visual_split=True),
    'UVVIS': SpectrumTechnique('UVVIS', x_reversed=False, threshold=0.05,
                          absorbance_label=True, visual_split=True,
                          em_wave=True),

    'THERMOGRAVIMETRIC ANALYSIS': SpectrumTechnique(
        'THERMOGRAVIMETRIC ANALYSIS', x_reversed=False, threshold=1.05),
    'X-RAY DIFFRACTION': SpectrumTechnique(
        'X-RAY DIFFRACTION', x_axis='xrd', x_reversed=False, threshold=1.00),
    'CYCLIC VOLTAMMETRY': SpectrumTechnique(
        'CYCLIC VOLTAMMETRY', x_axis='raw', y_axis='raw',
        x_reversed=False, threshold=1.00, cyclic_voltammetry=True),
    'SIZE EXCLUSION CHROMATOGRAPHY': SpectrumTechnique(
        'SIZE EXCLUSION CHROMATOGRAPHY', x_reversed=False, threshold=0.5,
        info_box='sec', sec_headers=True),
    'CIRCULAR DICHROISM SPECTROSCOPY': SpectrumTechnique(
        'CIRCULAR DICHROISM SPECTROSCOPY', x_reversed=False, threshold=1.00,
        negative_peaks=False),
    'SORPTION-DESORPTION MEASUREMENT': SpectrumTechnique(
        'SORPTION-DESORPTION MEASUREMENT', x_reversed=False, threshold=1.00),
    'Emissions': SpectrumTechnique('Emissions', x_reversed=False, threshold=0.5),
    'DLS ACF': SpectrumTechnique('DLS ACF', x_reversed=False, threshold=1.05),
    'DLS intensity': SpectrumTechnique('DLS intensity', x_reversed=False,
                                  threshold=1.00),

    # Forward, like TGA. It was reversed until #292, because `is_dsc` was
    # never added to the hand-maintained orientation chain in tf_img.
    'DIFFERENTIAL SCANNING CALORIMETRY': SpectrumTechnique(
        'DIFFERENTIAL SCANNING CALORIMETRY', x_reversed=False, threshold=1.05,
        info_box='dsc', dsc_metadata=True),

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

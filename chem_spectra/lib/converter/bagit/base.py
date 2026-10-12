import os
import base64
import tempfile
import json
import logging
import math

from chem_spectra.lib.converter.jcamp.base import (
    JcampBaseConverter, header_records, read_header,
)
from chem_spectra.lib.converter.share import (
    UnconvertibleSpectrum, parse_params,
)
from chem_spectra.lib.shared.misc import shorten_label
from chem_spectra.lib.converter.jcamp.data_parse import UnparsableJcampData
from chem_spectra.lib.converter.jcamp.technique import JcampTechniqueConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.composer.technique import (
    TechniqueComposer, flip_overlay_if_inverted,
)
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.lib.composer.lcms_converter_app import LCMSConverterAppComposer
from chem_spectra.lib.converter.share import parse_params
from chem_spectra.lib.converter.bagit.lcms_builder import append_lcms_group
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import ticker  # noqa: E402

logger = logging.getLogger(__name__)

# ##DATA TYPE spellings only an LC/MS run carries: its total ion chromatogram,
# or the LC/MS label ChemSpectra writes on the peak file it generates. Matched
# on the file's own text rather than on the mapped typ, because a request's
# data_type_mapping replaces the built-in mapping and need not list them.
LCMS_MARKER_DATATYPES = frozenset((
    'LC/MS', 'LCMS', 'LC-MS', 'MASS TIC',
    'TOTAL ION CHROMATOGRAM', 'TOTAL ION CHROMATOGRAPHY',
))
MASS_SPECTRUM_DATATYPES = frozenset(('MASS SPECTRUM', 'CONTINUOUS MASS SPECTRUM'))
HPLC_UVVIS_DATATYPES = frozenset(('HPLC UV/VIS SPECTRUM', 'HPLC UV-VIS'))
UVVIS_DATATYPES = frozenset(('UV/VIS SPECTRUM', 'UV-VIS', 'ULTRAVIOLET SPECTRUM'))


def _is_lcms_marker(cv):
    return cv.typ == 'LC/MS' or bool(set(cv.datatypes) & LCMS_MARKER_DATATYPES)


def _is_mass_spectrum(cv):
    return cv.typ == 'MS' or bool(set(cv.datatypes) & MASS_SPECTRUM_DATATYPES)


def _is_paged_ntuples(path):
    """Whether the block is an NTUPLES table read page by page.

    From the header text (#314's reader), since the base converter's
    dataclass only knows XYPOINTS and XYDATA.
    """
    records = list(header_records(read_header(path)))
    return (any(label == 'DATACLASS' and 'NTUPLES' in value.upper()
                for label, value in records)
            and any(label == 'PAGE' for label, _ in records))


def _is_lc_uv(cv, path):
    """Whether the member is a UV/VIS chromatogram, the LC half of a run.

    An HPLC UV/VIS block under either spelling, or a UV/VIS block that comes
    as wavelength pages. The converter's DATA TYPE is a profile dropdown, so
    the same LC-UV run can be labelled UV-VIS; its pages are what tells it
    from a spectrum, and only the LC/MS composer reads them.
    """
    datatypes = set(cv.datatypes)
    if cv.typ == 'HPLC UVVIS' or datatypes & HPLC_UVVIS_DATATYPES:
        return True
    is_uvvis = cv.typ == 'UVVIS' or bool(datatypes & UVVIS_DATATYPES)
    return is_uvvis and _is_paged_ntuples(path)


def _has_lcms_evidence(detected):
    """Whether the archive members together make an LC/MS run.

    A member is LC/MS outright (a TIC, or a re-uploaded LC/MS peak file), or
    is a UV/VIS chromatogram, or mass spectra come with a UV/VIS spectrum. A
    plain UV/VIS spectrum alone is not enough. `detected` maps each member's
    path to its base converter.
    """
    if any(_is_lcms_marker(cv) or _is_lc_uv(cv, path)
           for path, cv in detected.items()):
        return True
    converters = detected.values()
    has_ms = any(_is_mass_spectrum(cv) for cv in converters)
    return has_ms and any(cv.typ == 'UVVIS' for cv in converters)


class BagItBaseConverter:
    def __init__(self, target_dir, params=False, fname=''):
        self.raw_params = params
        self.params = parse_params(params)
        self.archive_entry_stems = []
        if target_dir is None:
            self.data, self.images, self.list_csv, self.combined_image = None, None, None, None
        else:
            ret = self.__read(target_dir, fname)
            if ret is None:
                self.data, self.images, self.list_csv, self.combined_image = None, None, None, None
            else:
                self.data, self.images, self.list_csv, self.combined_image = ret

    def __read(self, target_dir, fname):
        list_file_names = []
        data_dir_path = os.path.join(target_dir, 'data')
        flat_layout = not os.path.isdir(data_dir_path)
        if flat_layout:
            data_dir_path = target_dir
        for (dirpath, dirnames, filenames) in os.walk(data_dir_path):
            filenames.sort()
            list_file_names.extend(filenames)
            break
        if flat_layout:
            list_file_names = [n for n in list_file_names if n.lower().endswith('.jdx')]
        if (len(list_file_names) == 0):
            return None

        list_files = []
        list_images = []
        list_csv = []
        list_composer = []
        lcms_paths = []
        archive_stems = []
        # Every member is read once here and kept for the second pass, so
        # the archive is judged as a whole before anything is grouped.
        detected = {}
        for file_name in list_file_names:
            if not file_name.lower().endswith('.jdx'):
                continue
            jcamp_path = os.path.join(data_dir_path, file_name)
            try:
                detected[jcamp_path] = JcampBaseConverter(
                    jcamp_path, self.raw_params)
            except UnconvertibleSpectrum:
                # A member the converter refuses (a 2D file, say) refuses the
                # archive whole. Here, before any member is converted or drawn,
                # rather than from the loop below once earlier members have
                # already been rendered onto the shared figure.
                raise
            except Exception:
                pass
        has_lcms_context = _has_lcms_evidence(detected)

        for file_name in list_file_names:
            if not file_name.lower().endswith('.jdx'):
                continue
            jcamp_path = os.path.join(data_dir_path, file_name)
            stem = os.path.splitext(file_name)[0].replace('.', '_')
            base_cv = detected.get(jcamp_path) or JcampBaseConverter(
                jcamp_path, self.raw_params)
            # BagIt / flat LCMS zips: keep all chromatogram and MS traces in one
            # LCMSConverterAppComposer (incl. MASS SPECTRUM), not JcampMSConverter/ms.py.
            # Only an archive that is an LC/MS run: a UV/VIS spectrum on its
            # own (the converter ships every table as a BagIt) is a UV/VIS
            # spectrum, as it is when it arrives as a single file.
            # Also by the file's own DATA TYPE: a mass spectrum or chromatogram
            # the request's mapping does not name gets typ '' and would
            # otherwise go to the technique converter, which cannot read it.
            is_lcms_candidate = has_lcms_context and (
                base_cv.typ in ('HPLC UVVIS', 'UVVIS')
                or _is_lcms_marker(base_cv)
                or _is_mass_spectrum(base_cv)
                or _is_lc_uv(base_cv, jcamp_path))
            if is_lcms_candidate:
                lcms_paths.append(jcamp_path)
            else:
                try:
                    if base_cv.typ == 'MS':
                        # Standalone MS path
                        mscv = JcampMSConverter(base_cv)
                        tcp = MSComposer(mscv)
                    else:
                        try:
                            tcv = JcampTechniqueConverter(base_cv)
                        except UnparsableJcampData:
                            # one unusable member should not fail the whole
                            # archive; the rest still convert
                            logger.warning(
                                'no parsable data in %r inside the archive; '
                                'skipping it', jcamp_path,
                            )
                            continue
                        tcp = TechniqueComposer(tcv)
                except KeyError as err:
                    print(f"Skip empty JCAMP {file_name}: {err}")
                    continue
                # Carried on the composer rather than in a parallel list:
                # __combine_images filters LC/MS and MS composers out again, so
                # any index-aligned list would silently mislabel the survivors.
                tcp.source_filename = file_name
                list_composer.append(tcp)
                tf_jcamp = tcp.tf_jcamp()
                list_files.append(tf_jcamp)
                tf_img = tcp.tf_img()
                list_images.append(tf_img)
                if base_cv.typ == 'MS':
                    list_csv.append(None)
                else:
                    tf_csv = tcp.tf_csv()
                    list_csv.append(tf_csv)
                archive_stems.append(stem)

        if lcms_paths and parse_params(self.raw_params).get('transmittance'):
            # The LC/MS composer has no conversion, and these members are
            # read as one LC/MS dataset rather than as separate spectra, so
            # the instruction cannot be honoured here -- not even for a
            # UV/VIS member that converts perfectly well on its own. Said
            # rather than dropped: the archive came back 200 and unconverted,
            # which is indistinguishable from a conversion that happened.
            raise UnconvertibleSpectrum(
                'this archive is read as one LC/MS dataset, which has no '
                'transmittance conversion; convert the absorbance members '
                'on their own instead'
            )

        append_lcms_group(
            lcms_paths, self.raw_params,
            list_files, list_images, list_csv, list_composer,
            archive_stems=archive_stems,
        )

        self.archive_entry_stems = archive_stems
        self._composers = list_composer

        combined_image = self.__combine_images(list_composer)

        return list_files, list_images, list_csv, combined_image

    def get_base64_data(self):
        if self.data is None:
            return None
        list_jcamps = []
        for tf_jcamp in self.data:
            jcamp = base64.b64encode(tf_jcamp.read()).decode("utf-8")
            list_jcamps.append(jcamp)
        return list_jcamps

    @property
    def spc_type(self):
        """Derive the spectrum type label from the contained composers.

        Returns the actual type (e.g. ``'NMR SPECTRUM'``, ``'INFRARED
        SPECTRUM'``, ``'lcms'``) rather than a blanket ``'bagit'``, so
        callers can distinguish NMR / CV / UVVIS / LC-MS BagIt payloads.
        """
        composers = getattr(self, '_composers', None) or []
        if not composers:
            return 'bagit'

        types = set()
        for c in composers:
            if isinstance(c, LCMSConverterAppComposer):
                types.add('lcms')
            elif hasattr(c, 'core'):
                core = c.core
                if getattr(core, 'typ', None) == 'NMR':
                    types.add(getattr(core, 'ncl', 'NMR'))
                else:
                    types.add(getattr(core, 'typ', ''))

        if len(types) == 1:
            return types.pop()
        return 'bagit'

    def __combine_images(self, list_composer):
        non_lcms_techniques = [
            c for c in list_composer
            if not isinstance(c, LCMSConverterAppComposer)
            and not isinstance(c.core, JcampMSConverter)
        ]
        if len(non_lcms_techniques) <= 1:
            return None
        list_composer = non_lcms_techniques

        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        
        cv_mode = False
        cv_abs_max = 0.0
        for idx, composer in enumerate(list_composer):
            # `list_file_names` used to be a parameter here and was never
            # passed, so every BagIt legend read 0, 1, 2 ...
            filename = shorten_label(getattr(composer, 'source_filename', None)) \
                or str(idx)
            
            xs, ys = composer.core.xs, composer.core.ys
            y_values = ys
            if composer._technique().cyclic_voltammetry:
                cv_state = (
                    composer.core.params.get('cyclicvoltaSt')
                    or composer.core.params.get('cyclicvolta')
                    or composer.core.params.get('cyclic_volta')
                ) or {}
                if isinstance(cv_state, str):
                    try:
                        cv_state = json.loads(cv_state)
                    except Exception:
                        cv_state = {}
                cv_display = cv_state.get('cvDisplay') or {}
                if isinstance(cv_display, str):
                    try:
                        cv_display = json.loads(cv_display)
                    except Exception:
                        cv_display = {}
                try:
                    scale = float(cv_display.get('yScaleFactor', 1.0))
                except Exception:
                    scale = 1.0
                if scale != 1.0:
                    y_values = ys * scale
                cv_mode = True
                try:
                    cv_abs_max = max(cv_abs_max, float(np.max(np.abs(y_values))))
                except Exception:
                    pass
            marker = ''
            if composer._technique().sorption_branches:
                first_x, last_x = xs[0], xs[len(xs)-1]
                if first_x <= last_x:
                    filename = 'ADSORPTION'
                    marker = '^'
                else:
                    filename = 'DESORPTION'
                    marker = 'v'

            plt.plot(xs, y_values, label=filename, marker=marker)
            # PLOT label
            if composer._technique().x_axis == 'xrd':
                waveLength = composer.core.params['waveLength']
                label = "X ({}), WL={} nm".format(composer.core.label['x'], waveLength['value'], waveLength['unit'])    # noqa: E501
                plt.xlabel((label), fontsize=18)
            elif (composer._technique().cyclic_voltammetry):
                plt.xlabel("{}".format(composer.core.label['x']), fontsize=18)
            else:
                plt.xlabel("X ({})".format(composer.core.label['x']), fontsize=18)

            if (composer._technique().cyclic_voltammetry):
                plt.ylabel("{}".format(composer.core.label['y']), fontsize=18)
            else:
                plt.ylabel("Y ({})".format(composer.core.label['y']), fontsize=18)

        if flip_overlay_if_inverted(plt.gca(), [
            bool(getattr(c.core, 'draw_y_inverted', False))
            for c in list_composer
        ]):
            logger.info(
                'archive overlay mixes inverted and upright spectra; drawing '
                'it upright',
            )

        if cv_mode and cv_abs_max > 0:
            exp = int(math.floor(math.log10(cv_abs_max))) if cv_abs_max > 0 else 0
            base = (10.0 ** exp) if exp != 0 else 1.0
            ax = plt.gca()
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _:
                f"{(y / base):.3g}"
            ))
            ax.yaxis.get_offset_text().set_visible(False)
            if exp != 0:
                ax.text(
                    0.0, 1,
                    r"$\times 10^{%d}$" % exp,
                    transform=ax.transAxes,
                    ha='left', va='bottom',
                    fontsize=14,
                    clip_on=False
                )

        plt.legend()
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img

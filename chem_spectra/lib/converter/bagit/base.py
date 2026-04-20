import os
import base64
import tempfile
import json
import math

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.lib.composer.lcms_converter_app import LCMSConverterAppComposer
from chem_spectra.lib.converter.share import parse_params
from chem_spectra.lib.converter.bagit.lcms_builder import append_lcms_group
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import ticker  # noqa: E402


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
        for file_name in list_file_names:
            jcamp_path = os.path.join(data_dir_path, file_name)
            stem = os.path.splitext(file_name)[0].replace('.', '_')
            base_cv = JcampBaseConverter(jcamp_path, self.raw_params)
            # BagIt / flat LCMS zips: keep all chromatogram and MS traces in one
            # LCMSConverterAppComposer (incl. MASS SPECTRUM), not JcampMSConverter/ms.py.
            if base_cv.typ in ('MS', 'LC/MS', 'UVVIS', 'HPLC UVVIS'):
                lcms_paths.append(jcamp_path)
            else:
                try:
                    nicv = JcampNIConverter(base_cv)
                except KeyError as err:
                    print(f"Skip empty JCAMP {file_name}: {err}")
                    continue
                nicp = NIComposer(nicv)
                list_composer.append(nicp)
                tf_jcamp = nicp.tf_jcamp()
                list_files.append(tf_jcamp)
                tf_img = nicp.tf_img()
                list_images.append(tf_img)
                tf_csv = nicp.tf_csv()
                list_csv.append(tf_csv)
                archive_stems.append(stem)

        append_lcms_group(
            lcms_paths, self.raw_params,
            list_files, list_images, list_csv, list_composer,
            archive_stems=archive_stems,
        )

        self.archive_entry_stems = archive_stems

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

    def __combine_images(self, list_composer, list_file_names = None):
        if len(list_composer) <= 1:
            return None
        if any(isinstance(c, LCMSConverterAppComposer) for c in list_composer):
            return None
        if isinstance(list_composer[0].core, JcampMSConverter):
            return None

        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        
        cv_mode = False
        cv_abs_max = 0.0
        for idx, composer in enumerate(list_composer):
            filename = str(idx)
            if (list_file_names is not None) and idx < len(list_file_names):
                filename = list_file_names[idx]
            
            xs, ys = composer.core.xs, composer.core.ys
            y_values = ys
            if composer.core.is_cyclic_volta:
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
            if composer.core.is_aif:
                first_x, last_x = xs[0], xs[len(xs)-1]
                if first_x <= last_x:
                    filename = 'ADSORPTION'
                    marker = '^'
                else:
                    filename = 'DESORPTION'
                    marker = 'v'

            plt.plot(xs, y_values, label=filename, marker=marker)
            # PLOT label
            if (composer.core.is_xrd):
                waveLength = composer.core.params['waveLength']
                label = "X ({}), WL={} nm".format(composer.core.label['x'], waveLength['value'], waveLength['unit'])    # noqa: E501
                plt.xlabel((label), fontsize=18)
            elif (composer.core.is_cyclic_volta):
                plt.xlabel("{}".format(composer.core.label['x']), fontsize=18)
            else:
                plt.xlabel("X ({})".format(composer.core.label['x']), fontsize=18)

            if (composer.core.is_cyclic_volta):
                plt.ylabel("{}".format(composer.core.label['y']), fontsize=18)
            else:
                plt.ylabel("Y ({})".format(composer.core.label['y']), fontsize=18)
        
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

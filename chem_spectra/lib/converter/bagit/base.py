import os
import base64
import tempfile

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.converter.jcamp.ms import JcampMSConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.lib.composer.ms import MSComposer
from chem_spectra.lib.converter.share import parse_params
import matplotlib.pyplot as plt  # noqa: E402
import json

class BagItBaseConverter:
    def __init__(self, target_dir, params=False, fname=''):
        self.params = parse_params(params)
        if target_dir is None:
            self.data, self.images, self.list_csv, self.combined_image = None, None, None, None
        else:
            self.data, self.images, self.list_csv, self.combined_image = self.__read(target_dir, fname)

    def __read_metadata(self, target_dir):
        metadata_dir_path = os.path.join(target_dir, 'metadata/')
        metadata_json_filename = None
        for (dirpath, dirnames, filenames) in os.walk(metadata_dir_path):
            if len(filenames) > 0: metadata_json_filename = filenames[0]
            break
        
        auto_metadata = {}
        if metadata_json_filename is not None:
            json_path = os.path.join(metadata_dir_path, metadata_json_filename)
            with open(json_path) as json_file:
                try:
                    metadata = json.loads(json_file.read())
                    tables = metadata["tables"]
                    for table in tables:
                        table_filename = table["fileName"].replace("data/", "")
                        table_header = table["header"]
                        metadata_values = {}
                        for table_key in table_header.keys():
                            metadata_values[table_key.upper()] = table_header[table_key]
                        auto_metadata[table_filename] = metadata_values
                except Exception as e:
                    pass
        return auto_metadata

    def __read(self, target_dir, fname):
        list_file_names = []
        data_dir_path = os.path.join(target_dir, 'data')
        for (dirpath, dirnames, filenames) in os.walk(data_dir_path):
            filenames.sort()
            list_file_names.extend(filenames)
            break
        if (len(list_file_names) == 0):
            return None

        list_files = []
        list_images = []
        list_csv = []
        list_composer = []
        auto_metadata = self.__read_metadata(target_dir)
        for file_name in list_file_names:
            jcamp_path = os.path.join(data_dir_path, file_name)
            base_cv = JcampBaseConverter(jcamp_path)
            metadata = None
            if file_name in auto_metadata:
                metadata = auto_metadata[file_name]

            if base_cv.typ == 'MS':
                mscv = JcampMSConverter(base_cv)
                mscp = MSComposer(mscv)
                list_composer.append(mscp)
                tf_jcamp = mscp.tf_jcamp()
                list_files.append(tf_jcamp)
                tf_img = mscp.tf_img()
                list_images.append(tf_img)
                tf_csv = mscp.tf_csv()
                list_csv.append(tf_csv)
            else:
                nicv = JcampNIConverter(base_cv)
                if nicv.auto_metadata == None: nicv.auto_metadata = metadata
                nicp = NIComposer(nicv)
                list_composer.append(nicp)
                tf_jcamp = nicp.tf_jcamp()
                list_files.append(tf_jcamp)
                tf_img = nicp.tf_img()
                list_images.append(tf_img)
                tf_csv = nicp.tf_csv()
                list_csv.append(tf_csv)
        
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
        if isinstance(list_composer[0].core, JcampMSConverter):
            return None

        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        
        for idx, composer in enumerate(list_composer):
            filename = str(idx)
            if (list_file_names is not None) and idx < len(list_file_names):
                filename = list_file_names[idx]
            
            xs, ys = composer.core.xs, composer.core.ys
            marker = ''
            if composer.core.is_aif:
                first_x, last_x = xs[0], xs[len(xs)-1]
                if first_x <= last_x:
                    filename = 'ADSORPTION'
                    marker = '^'
                else:
                    filename = 'DESORPTION'
                    marker = 'v'

            plt.plot(xs, ys, label=filename, marker=marker)
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
        
        plt.legend()
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img

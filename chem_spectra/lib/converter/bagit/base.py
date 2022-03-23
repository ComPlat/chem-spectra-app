import os
import base64

from chem_spectra.lib.converter.jcamp.base import JcampBaseConverter
from chem_spectra.lib.converter.jcamp.ni import JcampNIConverter
from chem_spectra.lib.composer.ni import NIComposer
from chem_spectra.lib.converter.share import parse_params


class BagItBaseConverter:
    def __init__(self, target_dir, params=False, fname=''):
        self.params = parse_params(params)
        self.data, self.images = self.__read(target_dir, fname)

    def __read(self, target_dir, fname):
        list_file_names = []
        data_dir_path = os.path.join(target_dir, 'data')
        for (dirpath, dirnames, filenames) in os.walk(data_dir_path):
            list_file_names.extend(filenames)
            break
        if (len(list_file_names) == 0):
            return None

        list_files = []
        list_images = []
        for file_name in list_file_names:
            jcamp_path = os.path.join(data_dir_path, file_name)
            base_cv = JcampBaseConverter(jcamp_path)
            nicv = JcampNIConverter(base_cv)
            nicp = NIComposer(nicv)
            tf_jcamp = nicp.tf_jcamp()
            list_files.append(tf_jcamp)
            tf_img = nicp.tf_img()
            list_images.append(tf_img)
        return list_files, list_images

    def get_base64_data(self):
        if self.data is None:
            return None
        list_jcamps = []
        for tf_jcamp in self.data:
            jcamp = base64.b64encode(tf_jcamp.read()).decode("utf-8")
            list_jcamps.append(jcamp)
        return list_jcamps

import io
from chem_spectra.model.transformer import TransformerModel as TraModel
from chem_spectra.controller.helper.file_container import FileContainer
from werkzeug.datastructures import FileStorage

target_dir = './tests/fixtures/'
source_dir = 'source/'
result_dir = 'result/'

file_inte_mpy_jdx = 'CHI-224_10.jdx'
file_result_inte_mpy_jdx = 'result_CHI-224_10.jdx'

params = {'peaks_str': '', 'select_x': '', 'ref_name': '', 'ref_value': '', 'scan': 0, 'thres': 0.0, 'mass': 78.046950192, 'molfile': None,
          'clear': False, 'predict': '{}', 'ext': 'jdx', 'integration': '', 'multiplicity': '', 'fname': 'CHI-224_10.jdx', 'simulatenmr': False, 'waveLength': ''}


def test_ni_composer():
    with open(target_dir + source_dir + file_inte_mpy_jdx, 'rb') as f:
        file = FileStorage(f)
        file = FileContainer(file)
        molfile = FileContainer(None)
        composer, _ = TraModel(file, molfile=molfile,
                               params=params).to_composer()
        assert composer.title == 'CHI-224_10'


# def test_tf_jcamp():
#     with open(target_dir + result_dir + 'result_CHI-224_10.jdx', 'rb') as f:
#         file_test = f.readlines()

#     with open(target_dir + source_dir + file_inte_mpy_jdx, 'rb') as f:
#         file = FileStorage(f)
#         file = FileContainer(file)
#         molfile = FileContainer(None)
#         composer = TraModel(file, molfile=molfile, params=params).to_composer()
#         tf_jcamp = composer.tf_jcamp()
#         jcamp_file_content = tf_jcamp.readlines()
#         assert jcamp_file_content == file_test

def test_integrals():
    test_integration = ['(7.5077135686230382916, 7.5755592725546643251, 1.9512717547792841621)\n(7.4369542455041628415, 7.5006376363111479932, 1)\n(7.3374749618252721461, 7.4082342849441475963, 1.9811011720186100238)', '\n']
    with open(target_dir + source_dir + file_inte_mpy_jdx, 'rb') as f:
        file = FileStorage(f)
        file = FileContainer(file)
        molfile = FileContainer(None)
        composer, _ = TraModel(file, molfile=molfile,
                               params=params).to_composer()
        integration = composer.gen_integration_info()
        assert integration == test_integration


def test_multiplicity():
    test_multi = ['(1, 7.5077135686230382916, 7.5755592725546643251, 7.5429235558567508946, 1.9512717547792841621, 4, m, A)\n(2, 7.4369542455041628415, 7.5006376363111479932, 7.4632254127310462266, 1, 2, m, B)\n(3, 7.3374749618252721461, 7.4082342849441475963, 7.3728951424095283684, 1.9811011720186100238, 4, m, C)', '\n']
    with open(target_dir + source_dir + file_inte_mpy_jdx, 'rb') as f:
        file = FileStorage(f)
        file = FileContainer(file)
        molfile = FileContainer(None)
        composer, _ = TraModel(file, molfile=molfile,
                               params=params).to_composer()
        multiplicity = composer.gen_mpy_integ_info()
        assert multiplicity == test_multi

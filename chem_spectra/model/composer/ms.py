import tempfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from chem_spectra.model.composer.base import BaseComposer


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_SPECTRUM_EDIT = '$$ === CHEMSPECTRA SPECTRUM EDIT ===\n'
TEXT_MS_DATA_TABLE = '##DATA TABLE= (XY..XY), PEAKS\n' # '##XYDATA= (X++(Y..Y))\n'
TEXT_ASSIGN_AUTO = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS AUTO ===\n'
TEXT_ASSIGN_EDIT = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS EDIT ===\n'
TEXT_PEAK_ASSIGN = '##PEAK ASSIGNMENTS=(XYA)\n'


class MsComposer(BaseComposer):
    def __init__(self, core):
        super().__init__(core)
        self.title = core.fname
        self.meta = self.__compose()


    def __gen_headers_spectrum_orig(self):
        return [
            '\n',
            TEXT_SPECTRUM_ORIG,
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}\n'.format('MASS SPECTRUM'),
            '##DATA CLASS= NTUPLES\n',
            '##ORIGIN=\n',
            '##OWNER=\n',
            '##SPECTROMETER/DATA SYSTEM=\n',
            '##.SPECTROMETER TYPE={}\n'.format('TRAP'),
            '##.INLET={}\n'.format('GC'),
            '##.IONIZATION MODE={}\n'.format('EI+'),
            '##$SCANAUTOTARGET={}\n'.format(self.core.scan_auto_target),
            '##$SCANEDITTARGET={}\n'.format(''),
            '##$SCANCOUNT={}\n'.format(len(self.core.datatables)),
            '##$THRESHOLD={}\n'.format(0.05),
        ]


    def __gen_ntuples_begin(self):
        return ['##NTUPLES={}\n'.format('MASS SPECTRUM')]


    def __gen_ntuples_end(self):
        return ['##END NTUPLES={}\n'.format('MASS SPECTRUM')]


    def __gen_config(self):
        return [
            '##VAR_NAME= MASS, INTENSITY, RETENTION TIME\n',
            '##SYMBOL= X, Y, T\n',
            '##VAR_TYPE= INDEPENDENT, DEPENDENT, INDEPENDENT\n',
            '##VAR_FORM= AFFN, AFFN, AFFN\n',
            '##VAR_DIM= , , 3\n',
            '##UNITS= M/Z, RELATIVE ABUNDANCE, SECONDS\n',
            '##FIRST= , , 1\n',
            '##LAST= , , {}\n'.format(len(self.core.datatables)),
        ]


    def __gen_ms_spectra(self):
        msspcs = []
        for idx, dt in enumerate(self.core.datatables):
            msspc = [
                '##PAGE={}\n'.format(idx + 1),
                '##NPOINTS={}\n'.format(dt['pts']),
                TEXT_MS_DATA_TABLE,
            ]
            msspcs = msspcs + msspc + dt['dt']
        return msspcs


    def __compose(self):
        meta = []
        meta.extend(self.__gen_headers_spectrum_orig())

        meta.extend(self.__gen_ntuples_begin())
        meta.extend(self.__gen_config())
        meta.extend(self.__gen_ms_spectra())
        meta.extend(self.__gen_ntuples_end())

        meta.extend(self.gen_ending())
        return meta


    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        # PLOT data
        idx = self.core.scan_auto_target - 1
        spc = self.core.spectra[idx]
        plt.bar(spc[:, 0], spc[:, 1])

        # PLOT label
        plt.xlabel('X (m/z)', fontsize=18)
        plt.ylabel('Y (Relative Abundance)', fontsize=18)
        plt.grid(False)

        # Save
        tf = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf, format='png')
        tf.seek(0)
        plt.clf()
        plt.cla()
        return tf

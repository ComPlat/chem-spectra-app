import matplotlib
matplotlib.use('Agg')

import tempfile  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from chem_spectra.lib.composer.base import BaseComposer  # noqa: E402


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_MS_DATA_TABLE = '##DATA TABLE= (XY..XY), PEAKS\n'  # '##XYDATA= (X++(Y..Y))\n'  # noqa


class MSComposer(BaseComposer):
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
            '##.SPECTROMETER TYPE={}\n'.format(self.core.dic.get('SPECTROMETER TYPE', '')),  # TRAP     # noqa: E501
            '##.INLET={}\n'.format(self.core.dic.get('INLET', '')),  # GC
            '##.IONIZATION MODE={}\n'.format(self.core.dic.get('IONIZATION MODE', '')),  # EI+  # noqa: E501
            '##$CSCATEGORY=SPECTRUM\n',
            '##$CSSCANAUTOTARGET={}\n'.format(self.core.auto_scan),
            '##$CSSCANEDITTARGET={}\n'.format(
                self.core.edit_scan or self.core.auto_scan
            ),
            '##$CSSCANCOUNT={}\n'.format(len(self.core.datatables)),
            '##$CSTHRESHOLD={}\n'.format(self.core.thres / 100),
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

    def __prism(self, spc):
        blues_x, blues_y, greys_x, greys_y = [], [], [], []
        thres = 0
        if spc.shape[0] > 0:  # RESOLVE_VSMBNAN2
            thres = spc[:, 1].max() * (self.core.thres / 100)

        for pt in spc:
            x, y = pt[0], pt[1]
            if y >= thres:
                blues_x.append(x)
                blues_y.append(y)
            else:
                greys_x.append(x)
                greys_y.append(y)
        return blues_x, blues_y, greys_x, greys_y

    def prism_peaks(self):
        idx = (self.core.edit_scan or self.core.auto_scan) - 1
        spc = self.core.spectra[idx]
        return self.__prism(spc) + tuple([idx+1])

    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        # PLOT data
        blues_x, blues_y, greys_x, greys_y, _ = self.prism_peaks()
        plt.bar(greys_x, greys_y, width=0, edgecolor='#dddddd')
        plt.bar(blues_x, blues_y, width=0, edgecolor='#1f77b4')

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
    
    def tf_csv(self):
        return None

    def generate_nmrium(self):
        return None
    
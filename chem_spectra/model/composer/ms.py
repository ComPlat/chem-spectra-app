import tempfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from chem_spectra.model.composer.base import BaseComposer


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_SPECTRUM_EDIT = '$$ === CHEMSPECTRA SPECTRUM EDIT ===\n'
TEXT_DATA_TABLE = '##XYDATA= (XY..XY)\n' # '##XYDATA= (X++(Y..Y))\n'
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
            '##DATA CLASS=XYDATA\n',
            '##ORIGIN={}\n'.format(''), # TBD
            '##OWNER={}\n'.format(''), # TBD
            '##XUNITS={}\n'.format('M/Z'),
            '##YUNITS={}\n'.format('RELATIVE ABUNDANCE'),
            '##XFACTOR={}\n'.format('1'),
            '##YFACTOR={}\n'.format('1'),
            '##FIRSTX={}\n'.format(self.core.xs[0]),
            '##LASTX={}\n'.format(self.core.xs[-1]),
            '##MAXX={}\n'.format(self.core.xs.max()),
            '##MAXY={}\n'.format(self.core.ys.max()),
            '##MINX={}\n'.format(self.core.xs.min()),
            '##MINY={}\n'.format(self.core.ys.min())
        ]


    def __compose(self):
        meta = []
        meta.extend(self.gen_headers_root())

        meta.extend(self.__gen_headers_spectrum_orig())
        meta.extend(self.gen_spectrum_orig())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_ending())
        return meta


    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        # PLOT data
        plt.bar(self.core.xs, self.core.ys)

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

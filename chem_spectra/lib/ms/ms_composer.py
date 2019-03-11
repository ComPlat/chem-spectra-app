import tempfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_SPECTRUM_EDIT = '$$ === CHEMSPECTRA SPECTRUM EDIT ===\n'
TEXT_DATA_TABLE = '##XYDATA= (XY..XY)\n' # '##XYDATA= (X++(Y..Y))\n'
TEXT_ASSIGN_AUTO = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS AUTO ===\n'
TEXT_ASSIGN_EDIT = '$$ === CHEMSPECTRA PEAK ASSIGNMENTS EDIT ===\n'
TEXT_PEAK_ASSIGN = '##PEAK ASSIGNMENTS=(XYA)\n'


class MsComposer():
    def __init__(self, raw_converter):
        self.fname = raw_converter.fname
        self.xs, self.ys = self.__extract_spectrum(raw_converter)
        self.meta = self.__compose_block_meta()


    def __extract_spectrum(self, rc):
        spc = rc.spectrum
        return spc[:, 0], spc[:, 1]


    def __gen_headers_root(self):
        return [
            '##TITLE={}\n'.format(self.fname),
            '##JCAMP-DX=5.0\n',
            '##DATA TYPE=LINK\n',
            '##BLOCKS=1\n', # TBD
            '\n'
        ]


    def __header_base(self):
        return [
            '\n',
            TEXT_SPECTRUM_ORIG,
            '##TITLE={}\n'.format(self.fname),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}\n'.format('MASS SPECTRUM'),
            '##DATA CLASS=XYDATA\n',
            '##ORIGIN={}\n'.format(''), # TBD
            '##OWNER={}\n'.format(''), # TBD
        ]


    def __header_params(self):
        return [
            '##XUNITS={}\n'.format('M/Z'),
            '##YUNITS={}\n'.format('RELATIVE ABUNDANCE'),
            '##XFACTOR={}\n'.format('1'),
            '##YFACTOR={}\n'.format('1'),
            '##FIRSTX={}\n'.format(self.xs[0]),
            '##LASTX={}\n'.format(self.xs[-1]),
            '##MAXX={}\n'.format(self.xs.max()),
            '##MAXY={}\n'.format(self.ys.max()),
            '##MINX={}\n'.format(self.xs.min()),
            '##MINY={}\n'.format(self.ys.min())
        ]


    def __gen_headers_spectrum_orig(self):
        return self.__header_base() + self.__header_params()


    def __gen_datatable(self):
        pts = self.xs.shape[0]
        datatable = []
        for idx in range(pts):
            datatable.append(
                '{}, {}\n'.format(
                    self.xs[idx],
                    self.ys[idx]
                )
            )
        return datatable


    def __gen_spectrum_orig(self):
        c_spectrum_orig = [
            '##NPOINTS={}\n'.format(self.xs.shape[0]),
            TEXT_DATA_TABLE
        ]
        c_spectrum_orig.extend(self.__gen_datatable())
        return c_spectrum_orig


    def __gen_ending(self):
        return [
            '##END=\n',
            '\n'
        ]


    def __compose_block_meta(self):
        meta = []
        meta.extend(self.__gen_headers_root())

        meta.extend(self.__gen_headers_spectrum_orig())
        meta.extend(self.__gen_spectrum_orig())
        meta.extend(self.__gen_ending())
        meta.extend(self.__gen_ending())
        return meta


    def tf_jcamp(self):
        meta = ''.join(self.meta)
        tf = tempfile.NamedTemporaryFile(suffix='.jdx')
        with open(tf.name, 'w') as f:
            tf.write(bytes(meta, 'UTF-8'))
            tf.seek(0)
        return tf


    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        # PLOT data
        plt.bar(self.xs, self.ys)

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

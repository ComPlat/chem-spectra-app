import tempfile  # noqa: E402

from chem_spectra.lib.composer.base import BaseComposer  # noqa: E402
import numpy as np  # noqa: E402


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_MS_DATA_TABLE = '##DATA TABLE= (XY..XY), PEAKS\n'  # '##XYDATA= (X++(Y..Y))\n'  # noqa

class LCMSComposer:
    def __init__(self, core):
        self.core = core
        self.title = core.fname
        self.data = self.__compose()

    def __compose(self):
        tic_postive_data, tic_negative_data, uvvis_data, spectra_postive_data, spectra_negative_data = self.core.data
        tic_postive = self.__gen_tic(tic_data=tic_postive_data)
        tic_negative = self.__gen_tic(tic_data=tic_negative_data, is_negative=True)
        tic_positive_jcamp, tic_negative_jcamp = self.tf_jcamp(tic_postive), self.tf_jcamp(tic_negative)

        uvvis = self.__gen_uvvis(data=uvvis_data)
        uvvis_jcamp = self.tf_jcamp(uvvis)

        mz_positive = self.__gen_mz_spectra(data=spectra_postive_data)
        mz_positive_jcamp = self.tf_jcamp(mz_positive)

        mz_negative = self.__gen_mz_spectra(data=spectra_negative_data, is_negative=True)
        mz_negative_jcamp = self.tf_jcamp(mz_negative)
        return [tic_positive_jcamp, tic_negative_jcamp, uvvis_jcamp, mz_positive_jcamp, mz_negative_jcamp]

    def __gen_tic(self, tic_data, is_negative = False):
        time, intensity = tic_data['time'], tic_data['Intensity']
        max_time, min_time = np.max(time), np.min(time)
        max_intensity, min_intensity = np.max(intensity), np.min(intensity)
        content = [
            '\n',
            TEXT_SPECTRUM_ORIG,
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}\n'.format('LC/MS'),
            '##DATA CLASS= XYDATA\n',
            '##ORIGIN=\n',
            '##OWNER=\n',
            '##SPECTROMETER/DATA SYSTEM=\n',
            '##$CSCATEGORY=TIC SPECTRUM\n',
            '##XUNITS=time\n',
            '##YUNITS=Intensity\n',
            '##XFACTOR=1\n',
            '##YFACTOR=1\n',
            '##FIRSTX={}\n'.format(time[0]),
            '##LASTX={}\n'.format(time[len(time)-1]),
            '##MAXX={}\n'.format(max_time),
            '##MAXY={}\n'.format(max_intensity),
            '##MINX={}\n'.format(min_time),
            '##MINY={}\n'.format(min_intensity),
            '##NPOINTS={}\n'.format(len(time)),
            '##XYDATA= (XY..XY)\n',
        ]
        
        for i, _ in enumerate(time):
            content.append(
                '{}, {}\n'.format(time[i], intensity[i])
            )
        
        content.extend(self.__gen_ending())

        return content

    def __gen_uvvis(self, data):
        content = [
            '\n',
            TEXT_SPECTRUM_ORIG,
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}\n'.format('LC/MS'),
            '##DATA CLASS= NTUPLES\n',
            '##ORIGIN=\n',
            '##OWNER=\n',
            '##SPECTROMETER/DATA SYSTEM=\n',
            '##$CSCATEGORY=UVVIS SPECTRUM\n',
            '##VAR_NAME= RETENTION TIME, DETECTOR SIGNAL, WAVELENGTH\n',
            '##SYMBOL= X, Y, T\n',
            '##VAR_TYPE= INDEPENDENT, DEPENDENT, INDEPENDENT\n',
            '##VAR_FORM= AFFN, AFFN, AFFN\n',
            '##VAR_DIM= , , 3\n',
            '##UNITS= RETENTION TIME, DETECTOR SIGNAL, WAVELENGTH\n',
            '##FIRST= , , 1\n',
        ]
        
        msspcs = []
        ms_tempfile = tempfile.TemporaryFile()
        for time, value in data.items():
            xs, ys = value['RetentionTime'], value['DetectorSignal']
            msspc = [
                '##PAGE={}\n'.format(time),
                '##NPOINTS={}\n'.format(len(xs)),
                '##DATA TABLE= (XY..XY), PEAKS\n',
            ]
            for idx, _ in enumerate(xs):
                my_content = '{}, {};\n'.format(xs[idx], ys[idx])
                msspc += my_content
            file_content = ''.join(msspc)
            ms_tempfile.write(file_content.encode('utf-8'))

        ms_tempfile.seek(0)
        lines = ms_tempfile.readlines()
        decoded_lines = [line.decode('utf-8').strip() for line in lines]
        msspcs = '\n'.join(decoded_lines)
        ms_tempfile.close()
        
        content.extend(msspcs)
        content.extend(self.__gen_ending())

        return content

    def __gen_mz_spectra(self, data, is_negative=False):
        category = '##$CSCATEGORY=MZ NEGATIVE SPECTRUM\n' if is_negative else '##$CSCATEGORY=MZ POSITIVE SPECTRUM\n'
        content = [
            '\n',
            TEXT_SPECTRUM_ORIG,
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}\n'.format('LC/MS'),
            '##DATA CLASS= NTUPLES\n',
            '##ORIGIN=\n',
            '##OWNER=\n',
            '##SPECTROMETER/DATA SYSTEM=\n',
            '##$CSCATEGORY=UVVIS SPECTRUM\n',
            '##VAR_NAME= RETENTION TIME, DETECTOR SIGNAL, WAVELENGTH\n',
            '##SYMBOL= X, Y, T\n',
            '##VAR_TYPE= INDEPENDENT, DEPENDENT, INDEPENDENT\n',
            '##VAR_FORM= AFFN, AFFN, AFFN\n',
            '##VAR_DIM= , , 3\n',
            '##UNITS= RETENTION TIME, DETECTOR SIGNAL, WAVELENGTH\n',
            '##FIRST= , , 1\n',
            category,
        ]
        
        msspcs = []
        ms_tempfile = tempfile.TemporaryFile()
        for time, value in data.items():
            xs, ys = value['mz'], value['intensities']
            msspc = [
                '##PAGE={}\n'.format(time),
                '##NPOINTS={}\n'.format(len(xs)),
                '##DATA TABLE= (XY..XY), PEAKS\n',
            ]
            for idx, _ in enumerate(xs):
                my_content = '{}, {};\n'.format(xs[idx], ys[idx])
                msspc += my_content
            file_content = ''.join(msspc)
            ms_tempfile.write(file_content.encode('utf-8'))

        ms_tempfile.seek(0)
        lines = ms_tempfile.readlines()
        decoded_lines = [line.decode('utf-8').strip() for line in lines]
        msspcs = '\n'.join(decoded_lines)
        ms_tempfile.close()
        
        content.extend(msspcs)
        content.extend(self.__gen_ending())

        return content
        
    def __gen_ending(self):
        return [
            '##END=\n',
            '\n'
        ]
    
    def tf_jcamp(self, data):
        meta = ''.join(data)
        tf = tempfile.NamedTemporaryFile(suffix='.jdx')
        tf.write(bytes(meta, 'UTF-8'))
        tf.seek(0)
        return tf
# class LCMSComposer(BaseComposer):
#     def __init__(self, core):
#         super().__init__(core)
#         self.title = core.fname
#         self.meta = self.__compose()

    # def __gen_headers_spectrum_orig(self):
    #     return [
    #         '\n',
    #         TEXT_SPECTRUM_ORIG,
    #         '##TITLE={}\n'.format(self.title),
    #         '##JCAMP-DX=5.00\n',
    #         '##DATA TYPE={}\n'.format('LC/MS'),
    #         '##DATA CLASS= NTUPLES\n',
    #         '##ORIGIN=\n',
    #         '##OWNER=\n',
    #         '##SPECTROMETER/DATA SYSTEM=\n',
    #         # '##.SPECTROMETER TYPE={}\n'.format(self.core.dic.get('SPECTROMETER TYPE', '')),  # TRAP     # noqa: E501
    #         # '##.INLET={}\n'.format(self.core.dic.get('INLET', '')),  # GC
    #         # '##.IONIZATION MODE={}\n'.format(self.core.dic.get('IONIZATION MODE', '')),  # EI+  # noqa: E501
    #         '##$CSCATEGORY=SPECTRUM\n',
    #         # '##$CSSCANAUTOTARGET={}\n'.format(self.core.auto_scan),
    #         # '##$CSSCANEDITTARGET={}\n'.format(
    #         #     self.core.edit_scan or self.core.auto_scan
    #         # ),
    #         # '##$CSSCANCOUNT={}\n'.format(len(self.core.datatables)),
    #         # '##$CSTHRESHOLD={}\n'.format(self.core.thres / 100),
    #     ]

#     def __gen_ntuples_begin(self):
#         return ['##NTUPLES={}\n'.format('MASS SPECTRUM')]

#     def __gen_ntuples_end(self):
#         return ['##END NTUPLES={}\n'.format('MASS SPECTRUM')]

#     def __gen_config(self):
#         return [
            # '##VAR_NAME= MASS, INTENSITY, RETENTION TIME\n',
            # '##SYMBOL= X, Y, T\n',
            # '##VAR_TYPE= INDEPENDENT, DEPENDENT, INDEPENDENT\n',
            # '##VAR_FORM= AFFN, AFFN, AFFN\n',
            # '##VAR_DIM= , , 3\n',
            # '##UNITS= M/Z, RELATIVE ABUNDANCE, SECONDS\n',
            # '##FIRST= , , 1\n',
            # # '##LAST= , , {}\n'.format(len(self.core.datatables)),
#         ]

#     def __gen_ms_spectra(self):
        # msspcs = []
        # ms_tempfile = tempfile.TemporaryFile()
        # spectra_data = self.core.data[3] # the 1st and 2nd is tic positive and negative, the 3rd is uvvis
        # for time, value in spectra_data.items():
        #     xs, ys = value['mz'], value['intensities']
        #     msspc = [
        #         '##PAGE={}\n'.format(time),
        #         '##NPOINTS={}\n'.format(len(value['mz'])),
        #         '##DATA TABLE= (XY..XY), PEAKS\n',
        #     ]
        #     for idx in range(len(xs)):
        #         my_content = '{}, {};\n'.format(xs[idx], ys[idx])
        #         msspc += my_content
        #     file_content = ''.join(msspc)
        #     ms_tempfile.write(file_content.encode('utf-8'))

        # ms_tempfile.seek(0)
        # lines = ms_tempfile.readlines()
        # decoded_lines = [line.decode('utf-8').strip() for line in lines]
        # msspcs = '\n'.join(decoded_lines)
        # ms_tempfile.close()
#         return msspcs

#     def __compose(self):
#         meta = []
#         meta.extend(self.__gen_headers_spectrum_orig())

#         meta.extend(self.__gen_ntuples_begin())
#         meta.extend(self.__gen_config())
#         meta.extend(self.__gen_ms_spectra())
#         meta.extend(self.__gen_ntuples_end())

#         # meta.extend(self.generate_original_metadata())

#         meta.extend(self.gen_ending())
#         return meta

#     # def __prism(self, spc):
#     #     blues_x, blues_y, greys_x, greys_y = [], [], [], []
#     #     thres = 0
#     #     if spc.shape[0] > 0:  # RESOLVE_VSMBNAN2
#     #         thres = spc[:, 1].max() * (self.core.thres / 100)

#     #     for pt in spc:
#     #         x, y = pt[0], pt[1]
#     #         if y >= thres:
#     #             blues_x.append(x)
#     #             blues_y.append(y)
#     #         else:
#     #             greys_x.append(x)
#     #             greys_y.append(y)
#     #     return blues_x, blues_y, greys_x, greys_y

#     # def prism_peaks(self):
#     #     idx = (self.core.edit_scan or self.core.auto_scan) - 1
#     #     spc = self.core.spectra[idx]
#     #     return self.__prism(spc) + tuple([idx+1])

#     def tf_img(self):
#         # plt.rcParams['figure.figsize'] = [16, 9]
#         # plt.rcParams['font.size'] = 14
#         # # PLOT data
#         # blues_x, blues_y, greys_x, greys_y, _ = self.prism_peaks()
#         # plt.bar(greys_x, greys_y, width=0, edgecolor='#dddddd')
#         # plt.bar(blues_x, blues_y, width=0, edgecolor='#1f77b4')

#         # # PLOT label
#         # plt.xlabel('X (m/z)', fontsize=18)
#         # plt.ylabel('Y (Relative Abundance)', fontsize=18)
#         # plt.grid(False)

#         # # Save
#         # tf = tempfile.NamedTemporaryFile(suffix='.png')
#         # plt.savefig(tf, format='png')
#         # tf.seek(0)
#         # plt.clf()
#         # plt.cla()
#         # return tf
#         return None
    
#     def tf_csv(self):
#         return None

#     def generate_nmrium(self):
#         return None
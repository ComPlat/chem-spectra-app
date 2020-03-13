import matplotlib
matplotlib.use('Agg')

import tempfile  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.path as mpath  # noqa: E402
import re # noqa: E402
import numpy as np # noqa: E402

from chem_spectra.lib.composer.base import (  # noqa: E402
    extrac_dic, calc_npoints, BaseComposer
)
from chem_spectra.lib.shared.calc import (  # noqa: E402
    calc_mpy_center, calc_ks, get_curve_endpoint
)


TEXT_SPECTRUM_ORIG = '$$ === CHEMSPECTRA SPECTRUM ORIG ===\n'
TEXT_INTEGRATION = '$$ === CHEMSPECTRA INTEGRATION ===\n'
TEXT_MULTIPLICITY = '$$ === CHEMSPECTRA MULTIPLICITY ===\n'


class NIComposer(BaseComposer):
    def __init__(self, core):
        super().__init__(core)
        self.title = core.fname
        self.meta = self.__compose()

    def __header_base(self):
        return [
            '\n',
            TEXT_SPECTRUM_ORIG,
            '##TITLE={}\n'.format(self.title),
            '##JCAMP-DX=5.00\n',
            '##DATA TYPE={}\n'.format(self.core.datatype),
            '##DATA CLASS=XYDATA\n',
            '##$CSCATEGORY=SPECTRUM\n',
            '##ORIGIN={}\n'.format(extrac_dic(self.core, 'ORIGIN')),
            '##OWNER={}\n'.format(extrac_dic(self.core, 'OWNER')),
        ]

    def __calc_nucleus_by_boundary(self):
        delta = self.core.boundary['x']['max'] - self.core.boundary['x']['min']
        nucleus = '^1H' if abs(delta) < 40.0 else '^13C'
        return nucleus


    def __get_nucleus(self):
        nuc_orig = extrac_dic(self.core, '.OBSERVENUCLEUS')
        nuc_modf = re.sub('[^A-Za-z0-9]+', '', nuc_orig).lower()
        is_valid = ('13c' in nuc_modf) or ('1h' in nuc_modf) or ('19f' in nuc_modf)
        nucleus = nuc_orig if is_valid else self.__calc_nucleus_by_boundary()
        return nucleus

    def __header_nmr(self):
        return [
            '##.OBSERVE FREQUENCY={}\n'.format(
                extrac_dic(self.core, '.OBSERVEFREQUENCY')
            ),
            '##.OBSERVE NUCLEUS={}\n'.format(
                self.__get_nucleus()
            ),
            '##SPECTROMETER/DATA SYSTEM={}\n'.format(
                extrac_dic(self.core, 'SPECTROMETER/DATASYSTEM')
            ),
        ]

    def __header_params(self):
        return [
            '##XUNITS={}\n'.format(self.core.label['x']),
            '##YUNITS={}\n'.format(self.core.label['y']),
            '##XFACTOR={}\n'.format(self.core.factor['x']),
            '##YFACTOR={}\n'.format(self.core.factor['y']),
            '##FIRSTX={}\n'.format(self.core.boundary['x']['max']),
            '##LASTX={}\n'.format(self.core.boundary['x']['min']),
            '##MAXX={}\n'.format(self.core.boundary['x']['max']),
            '##MAXY={}\n'.format(self.core.boundary['y']['max']),
            '##MINX={}\n'.format(self.core.boundary['x']['min']),
            '##MINY={}\n'.format(self.core.boundary['y']['min'])
        ]

    def __gen_headers_spectrum_orig(self):
        if self.core.typ == 'INFRARED':
            return self.__header_base() + self.__header_params()
        else:
            return self.__header_base() + \
                self.__header_nmr() + self.__header_params()

    def __gen_headers_im(self):
        return [
            '$$ === CHEMSPECTRA ===\n',
        ]

    def __gen_headers_integration(self):
        return [
            '##$OBSERVEDINTEGRALS= (X Y Z)\n',
        ]

    def __gen_headers_mpy_integ(self):
        return [
            '##$OBSERVEDMULTIPLETS=\n',
        ]

    def __gen_headers_mpy_peaks(self):
        return [
            '##$OBSERVEDMULTIPLETSPEAKS=\n',
        ]

    def __compose(self):
        meta = []
        meta.extend(self.gen_headers_root())

        meta.extend(self.__gen_headers_spectrum_orig())
        meta.extend(self.gen_spectrum_orig())
        meta.extend(self.__gen_headers_im())
        meta.extend(self.__gen_headers_integration())
        meta.extend(self.gen_integration_info())
        meta.extend(self.__gen_headers_mpy_integ())
        meta.extend(self.gen_mpy_integ_info())
        meta.extend(self.__gen_headers_mpy_peaks())
        meta.extend(self.gen_mpy_peaks_info())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_headers_peaktable_edit())
        meta.extend(self.gen_edit_peaktable())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_headers_peaktable_auto())
        meta.extend(self.gen_auto_peaktable())
        meta.extend(self.gen_ending())

        meta.extend(self.gen_ending())
        return meta

    def __plt_nbins(self):
        typ = self.core.typ
        if 'NMR' == typ:
            return 20
        elif 'INFRARED' == typ:
            return 20
        elif 'MS' == typ:
            return 20
        return 20

    def __fakto(self):
        typ = self.core.typ
        if 'NMR' == typ:
            return 1
        elif 'MS' == typ:
            return 1
        elif 'INFRARED' == typ:
            return -1
        return 1

    def tf_img(self):
        plt.rcParams['figure.figsize'] = [16, 9]
        plt.rcParams['font.size'] = 14
        # PLOT data
        plt.plot(self.core.xs, self.core.ys)
        plt.xlim(
            self.core.boundary['x']['max'],
            self.core.boundary['x']['min']
        )
        y_max, y_min = np.max(self.core.ys), np.min(self.core.ys)
        h = y_max - y_min
        plt.ylim(
            y_min - h * 0.2,
            y_max + h * 0.2,
        )
        # PLOT peaks
        faktor = self.__fakto()
        path_data = [
            (mpath.Path.MOVETO, (0, faktor * 5)),
            (mpath.Path.LINETO, (0, faktor * 20)),
        ]
        codes, verts = zip(*path_data)
        marker = mpath.Path(verts, codes)
        if self.core.edit_peaks:
            plt.plot(
                self.core.edit_peaks['x'],
                self.core.edit_peaks['y'],
                'rd',
                marker=marker,
                markersize=50,
            )
        elif self.core.auto_peaks:
            plt.plot(
                self.core.auto_peaks['x'],
                self.core.auto_peaks['y'],
                'rd',
                marker=marker,
                markersize=50,
            )

        # ----- PLOT integration -----
        refShift, refArea = self.refShift, self.refArea
        itg_h = y_max + h * 0.1
        for itg in self.all_itgs:
            # integration marker
            xL, xU, area = itg['xL'] - refShift, itg['xU'] - refShift, itg['area'] * refArea
            plt.plot([xL, xU], [itg_h, itg_h], color='#228B22')
            plt.plot([xL, xL], [itg_h + h * 0.01, itg_h - h * 0.01], color='#228B22')
            plt.plot([xU, xU], [itg_h + h * 0.01, itg_h - h * 0.01], color='#228B22')
            plt.text((xL + xU) / 2, itg_h + h * 0.015, '{:0.2f}'.format(area), color='#228B22', ha='center', size=12)
            # integration curve
            ks = calc_ks(self.core.ys, y_max, h)
            iL, iU = get_curve_endpoint(self.core.xs, self.core.ys, xL, xU)
            ref = ks[iL]
            cxs = self.core.xs[iL:iU]
            cys = (ks[iL:iU] - ref) * 1.5 + (y_max - h * 0.4)
            plt.plot(cxs, cys, color='#228B22')

        # ----- PLOT multiplicity -----
        mpy_h = y_min - h * 0.08
        for mpy in self.mpys:
            xL, xU, area, typ, peaks = mpy['xExtent']['xL'] - refShift, mpy['xExtent']['xU'] - refShift, mpy['area'] * refArea, mpy['mpyType'], mpy['peaks']
            plt.plot([xL, xU], [mpy_h, mpy_h], color='#DA70D6')
            plt.plot([xL, xL], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')
            plt.plot([xU, xU], [mpy_h + h * 0.01, mpy_h - h * 0.01], color='#DA70D6')
            plt.text((xL + xU) / 2, mpy_h - h * 0.1, '({})'.format(typ), color='#DA70D6', ha='center', size=12)
            plt.text((xL + xU) / 2, mpy_h - h * 0.06, '{:0.3f}'.format(calc_mpy_center(mpy['peaks'], refShift, mpy['mpyType'])), color='#DA70D6', ha='center', size=12)
            for p in peaks:
                x = p['x']
                plt.plot([x - refShift, x - refShift], [mpy_h, mpy_h + h * 0.05], color='#DA70D6')

        # PLOT label
        plt.xlabel("X ({})".format(self.core.label['x']), fontsize=18)
        plt.ylabel("Y ({})".format(self.core.label['y']), fontsize=18)
        plt.locator_params(nbins=self.__plt_nbins())
        plt.grid(False)
        # Save
        tf_img = tempfile.NamedTemporaryFile(suffix='.png')
        plt.savefig(tf_img, format='png')
        tf_img.seek(0)
        plt.clf()
        plt.cla()
        return tf_img

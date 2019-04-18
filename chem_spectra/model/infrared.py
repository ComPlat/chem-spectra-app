import numpy as np
import tempfile
from scipy import interpolate as sc_interpolate

from chem_spectra.model.helper.share import store_str_in_tmp
from chem_spectra.model.converter.jcamp.base import JcampBaseConverter
from chem_spectra.model.converter.jcamp.ni import JcampNIConverter


class InfraredModel:
    def __init__(self, spectrum):
        self.nicv = self.__read_spectrum(spectrum)


    def __order(self, x_i, y_i):
        x_o = x_i
        y_o = y_i

        first_x = x_i[0]
        last_x = x_i[-1]
        if first_x > last_x:
            x_o = x_o[::-1]
            y_o = y_o[::-1]
        return x_o, y_o


    def __concat_boundary(self, x_i, y_i):
        x_head = np.array([0.0])
        x_tail = np.array([4000.0])
        x_o = np.concatenate([x_head, x_i, x_tail])

        y_head = np.array([y_i[0]])
        y_tail = np.array([y_i[-1]])
        y_o = np.concatenate([y_head, y_i, y_tail])
        return x_o, y_o


    def __interpolate(self, x_i, y_i):
        f = sc_interpolate.interp1d(x_i, y_i)
        x_o = np.linspace(0, 3999, 4000, endpoint=True)
        y_o = f(x_o)
        return x_o, y_o


    def __normalize(self, x_i, y_i):
        max_y = np.max(y_i)
        min_y = np.min(y_i)
        height = max_y - min_y

        x_o = x_i
        y_o = (y_i - min_y) / height
        return x_o, y_o


    def __use_absorption(self, x_i, y_i):
        x_o = x_i
        y_o = y_i

        y_median = np.median(y_i)
        if y_median < 0.5:
            y_o = 1 - y_i
        return x_o, y_o


    def __chop(self, x_i, y_i):
        x_o = x_i[0:]
        y_o = y_i[0:]

        return x_o, y_o


    def __read_spectrum(self, spc):
        tf = store_str_in_tmp(spc.core, suffix='.jdx')
        jbcv = JcampBaseConverter(tf.name)
        nicv = JcampNIConverter(jbcv)
        tf.close()
        return nicv


    def standarize(self):
        x_order, y_order = self.__order(self.nicv.xs, self.nicv.ys)
        x_full, y_full = self.__concat_boundary(x_order, y_order)
        x_itp, y_itp = self.__interpolate(x_full, y_full)
        x_itp_norm, y_itp_norm = self.__normalize(x_itp, y_itp)
        x_itp_norm_abs, y_itp_norm_abs = self.__use_absorption(x_itp_norm, y_itp_norm)
        x_inac, y_inac = self.__chop(x_itp_norm_abs, y_itp_norm_abs)

        return x_inac, y_inac

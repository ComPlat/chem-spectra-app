# Reference
# https://github.com/jjhelmus/nmrglue/blob/master/nmrglue/fileio/jcampdx.py#L125
DICT_DIGITS = {"0": "@",
               "1": "A", "2": "B", "3": "C", "4": "D", "5": "E",
               "6": "F", "7": "G", "8": "H", "9": "I",
               "-1": "a", "-2": "b", "-3": "c", "-4": "d", "-5": "e",
               "-6": "f", "-7": "g", "-8": "h", "-9": "i"}


class DatatableModel:
    def __to_encode_str(self, num):
        num_abs_lst = list(str(abs(int(num))))
        if num < 0:
            num_abs_lst[0] = DICT_DIGITS['-' + num_abs_lst[0]]
        else:
            num_abs_lst[0] = DICT_DIGITS[num_abs_lst[0]]
        return ''.join(num_abs_lst)

    def __encode_row(self, x, ys):
        output = [str(x)]
        for idx, y in enumerate(ys):
            output.append(self.__to_encode_str(y))

        return ''.join(output) + '\n'

    def encode(self, arr_ys, y_factor, arr_xs=None, is_xypoints=False):
        total_count = arr_ys.shape[0]
        calculate_ys = []

        if (is_xypoints):
            for idx in range(len(arr_xs)):
                line = '{x}, {y}\n'.format(x=arr_xs[idx], y=arr_ys[idx])
                calculate_ys.append(line)
        else:
            ys = []
            line = None
            line_x = None
            for i, v in enumerate(arr_ys):
                v = int(v / float(y_factor))
                ys.append(v)
                if line is None:
                    line = str(i)
                    line_x = i
                line += str(abs(v))
                if len(line) > 66 or i+1 == total_count:
                    output = self.__encode_row(total_count - line_x, ys)
                    calculate_ys.append(output)
                    ys = []
                    line = None
                    line_x = None

        return calculate_ys

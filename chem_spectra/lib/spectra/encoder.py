import math


# Reference
# https://github.com/jjhelmus/nmrglue/blob/master/nmrglue/fileio/jcampdx.py#L125
DICT_DIGITS = {"0": "@",
               "1": "A", "2": "B", "3": "C", "4": "D", "5": "E",
               "6": "F", "7": "G", "8": "H", "9": "I",
               "-1": "a", "-2": "b", "-3": "c", "-4": "d", "-5": "e",
               "-6": "f", "-7": "g", "-8": "h", "-9": "i"}


def abs_digit_count(target):
    return len(str(abs(int(target))))


def to_encode_str(num):
    num_abs_lst = list(str(abs(int(num))))
    if num < 0:
        num_abs_lst[0] = DICT_DIGITS['-' + num_abs_lst[0]]
    else:
        num_abs_lst[0] = DICT_DIGITS[num_abs_lst[0]]
    return ''.join(num_abs_lst)


def encode_row(x, ys):
    output = [str(x)]
    for idx, y in enumerate(ys):
        output.append(to_encode_str(y))

    return ''.join(output) + '\n'


def encode_datatable(arr_ys, max_y, y_factor):
    total_count = arr_ys.shape[0]
    calculate_ys = []

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
            output = encode_row(total_count - line_x, ys)
            calculate_ys.append(output)
            ys = []
            line = None
            line_x = None

    return calculate_ys

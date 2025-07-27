#  Copyright (c) 2025. By ZYF

import os
import chardet

def format_file(file, to_format ='unix2dos'):
    print('Formatting %s: \t%s' % (to_format, file))
    if not os.path.isfile(file):
        raise FileNotFoundError("Not a file.")

    if to_format == 'unix2dos':
        line_sep = '\r\n'
    elif to_format == 'dos2unix':
        line_sep = '\n'
    else:
        raise ValueError('Invalid Parameter "toformat"')

    with open(file, 'rb') as f:
        raw_data = f.read()
        result = chardet.detect(raw_data)
        encoding = result['encoding']
        print('Encoding:\t' + encoding)

    with open(file, 'r', encoding=encoding) as fd:
        tmp_file = open(file + to_format, 'w+b')
        for line in fd:
            line = line.replace('\r', '')
            line = line.replace('\n', '')
            tmp_file.write((line + line_sep).encode())
        tmp_file.close()
    os.remove(file)
    os.renames(file + to_format, file)

    print('Successful formatting!')

if __name__ == '__main__':
    file = 'D:/Data/大学学习/计算化学与人工智能/苏亚琼课题组工作/毕业设计FeN4X/VASP input/Eform/Br/INCAR'
    format_file(file, 'dos2unix')
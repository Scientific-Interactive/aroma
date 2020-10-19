#
# Author: V. Ganesh
# Date: 19.10.2020
# Python 3 replacement for fpformat

def fix(x, digits):
    return '%' + repr(digits) + 'f' % x

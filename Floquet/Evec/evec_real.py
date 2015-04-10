__author__ = 'liangshengzhang'

"""
Check whether eigenstates are real up to an overall phase.
"""

path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing_new_1 as Draw
import label_build as Label
import numpy as np
from math import sin, cos

include_s = "Shift" # A string that must be included in filename construct
exclude_s = None # A string that must be excluded in filename construct
#ver_num = 1

def read_file(filename):
    temp = np.loadtxt(filename+".txt", dtype=complex)
    return temp

filename = []
label = [] # General label for excuding when plotting

# Obtain filenames
Label.name_read("name", filename, label, include_s = include_s, exclude_s=exclude_s)

legend = ['' for n in filename] # labels used as legends in plotting

# Generic labels
general = [[] for n in filename]

Label.string_extract(filename, general, end_string = ",size=")
Label.label_change(general)

#print general

# Length
length = [[] for n in filename]
length_start = "size="
length_end = ","
length_label = "L="

Label.string_extract(filename, length, start_string = length_start, end_string = length_end)

# Number of points
J = [[] for n in filename]
J_string = "J="

Label.string_extract(filename, J, start_string = J_string, end_string = ',')

# Version number if exists
version = [[] for n in filename]
version_start = ",v"
version_label = "v"

Label.string_extract(filename, version, start_string = version_start)

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(legend, '', general, end = True, concat = '')

Label.label_build(label, length_label, length, end = True)
Label.label_build(legend, length_label, length, end = True)

Label.label_build(label, J_string, J, end = True)
Label.label_build(legend, J_string, J, end = True)

Label.label_build(label, version_label, version, end = True)

# Remove XXZ part from legends
remove_string = "XXZ "
for n in range(len(legend)):
    s = legend[n]
    start = s.find(remove_string)
    if start > 0:
        legend[n] = s[:start-1]
    else:
        legend[n] = ''
    legend[n] += s[start+len(remove_string):]

index = 0
evec = []
for n in filename:
    print "Read file", n
    evec.append( read_file(filename[index]) )
    index += 1

# Check norm
for i in range(len(evec)):
    print "Check normalization case", i
    for j in range(len(evec[i])):
        s = np.linalg.norm(evec[i][j])
        if abs(s-1) > 1.0e-6:
            print "ver: ",i, j, s

# Check whether they are real
for l in range(len(evec)):
    print "Check real evec case", l
    for i in range(len(evec)):
        index = 0
        while np.linalg.norm(evec[l][i][index]) < 1.0e-6:
            index += 1
        s = np.angle(evec[l][i][index])
        scale = complex(cos(s), -sin(s))
        temp_evec = [scale * n for n in evec[l][i]]

        is_real = True
        first_imag = 0
        for j in range(len(temp_evec)):
            if abs( np.imag(temp_evec[j]) ) > 2 * 1.0e-6:
                is_real = False
                first_imag = np.imag(temp_evec[j])
                break

        if not is_real:
            print i, "Not real", first_imag

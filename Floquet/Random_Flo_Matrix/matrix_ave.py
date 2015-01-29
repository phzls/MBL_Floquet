path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

"""
This file averages matrices elements along diagonal direction.
"""

import drawing as Draw
import label_build as Label
import numpy as np
from math import sqrt

def read_matrix_norm(filename):
    """ Read raw matrix elements"""
    matrix = np.loadtxt(filename+".txt")
    return matrix

filename = []
label = []

# Obtain filenames
Label.name_read("name", filename, label)

legend = ['' for n in filename] # labels used as legends in plotting

# Generic labels
general = [[] for n in filename]

Label.string_extract(filename, general, end_string = "_L=")
Label.label_change(general)

# Length
length = [[] for n in filename]
length_start = "L="
length_end = ","
length_label = "L="

Label.string_extract(filename, length, start_string = length_start, end_string = length_end)

# Type of calculation
cal_type = [[] for n in filename]
type_string = "right_most_sigma_z_"
type_label = "Type="

Label.string_extract(filename, cal_type, start_string = type_string)
Label.label_change(cal_type)

# Realizations
realization = [[] for n in filename]
realization_string = "Realizations="
realization_label = "Run="

Label.string_extract(filename, realization, start_string = realization_string, end_string = ',')

# minimum angle for rotation if exists
angle_min = [[] for n in filename]
angle_min_start = "angle_min="
angle_min_label = "angle min="

Label.string_extract(filename, angle_min, start_string = angle_min_start, end_string = ',')

# Supreme angle for rotation if exists
angle_sup = [[] for n in filename]
angle_sup_start = "angle_sup="
angle_sup_label = "angle="

Label.string_extract(filename, angle_sup, start_string = angle_sup_start, end_string = ',')

# Coupling strength if exists
coupling = [[] for n in filename]
coupling_start = "J="
coupling_label = "J="

Label.string_extract(filename, coupling, start_string = coupling_start, end_string = ',')

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(label, length_label, length, end = True)
Label.label_build(label, coupling_label, coupling, end = True)
#Label.label_build(label, angle_min_label, angle_min, end = True)
Label.label_build(label, angle_sup_label, angle_sup, end = True)
Label.label_build(label, realization_label, realization, end = True)
Label.label_build(label, type_label,cal_type, end = True)
print label

Label.label_build(legend, '', general, end = True, concat = '')
Label.label_build(legend, length_label, length, end = True)
Label.label_build(legend, coupling_label, coupling, end = True)
#Label.label_build(label, angle_min_label, angle_min, end = True)
Label.label_build(legend, angle_sup_label, angle_sup, end = True)

matrix = []

for i in range(len(filename)):
    print "Read file", filename[i]
    matrix.append([])

    if cal_type[i] == "norm":
        matrix[i] = read_matrix_norm(filename[i])

diag_ave = []
diff = []
ave = []
std = []

for i in range(len(matrix)):
    diag_ave.append({})
    diff.append([])
    ave.append([])
    std.append([])
    if len(matrix[i]) != 0:
        for row in range(len(matrix[i])):
            for col in range(len(matrix[i][row])):
                diff_val = abs(row - col)
                if diff_val > len(matrix[i][col])/2:
                    diff_val = len(matrix[i][col]) - diff_val
                try:
                    diag_ave[i][diff_val]
                except KeyError:
                    diag_ave[i][diff_val] = []
                diag_ave[i][diff_val].append(matrix[i][row][col])

        for key in sorted(diag_ave[i].keys()):
            diff[i].append(key)
            ave[i].append(np.mean(diag_ave[i][key]))
            if len(diag_ave[i][key]) == 1:
                std[i].append(0)
            else:
                std[i].append(np.std(diag_ave[i][key], ddof=1) / sqrt(len(diag_ave[i][key])))
    print i, len(diag_ave[i]), diff[i][1020:1030], std[i][1020:1030]

import pylab
draw1 = Draw.Draw()

draw1.figure_init(logy = True)
draw1.figure_set()

in_range = ["Random"] # Models
must_in_range = ["norm", "Run=100 "]
must_not_range = ["J=0.9", "J=0.7"]

extra_index= [] # The XXZ model

#for n in xrange(len(label)):
#    if label[n].find("norm") > -1 and label[n].find("XXZ") > -1:
#        extra_index.append(n)

print extra_index

draw1.plot_range(label, in_range = in_range, must_in_range = must_in_range,
                  must_not_range = must_not_range, index = extra_index, printout = True)

draw1.errorbar(diff, ave, yerr = std, label = legend)

pylab.legend(loc='upper right', ncol=1, prop={'size':15}, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$\langle\sigma^z_{i,j} \rangle$")


#pylab.savefig("High_Disorder_Floquet_10_rightmost_sigma_z_entry_norm_ave_diag_run_100.pdf",
#    box_inches='tight')

pylab.show()


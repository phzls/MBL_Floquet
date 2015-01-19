path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

"""
This file gives histogram of matrices elements without diagonal elements for one realization.
"""

import drawing as Draw
import label_build as Label
import numpy as np

def read_matrix_norm(filename):
    """ Read raw level statistitcs"""
    matrix = np.loadtxt(filename+".txt")
    return matrix


filename = []
label = []

# Obtain filenames
Label.name_read("name", filename, label)

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

Label.string_extract(filename, cal_type, start_string = type_string)
Label.label_change(cal_type)

# Realizations
realization = [[] for n in filename]
realization_string = "Realizations="
realization_label = "Run="

Label.string_extract(filename, realization, start_string = realization_string, end_string = ',')
print realization

# minimum angle for rotation if exists
angle_min = [[] for n in filename]
angle_min_start = "angle_min="
angle_min_label = "angle min="

Label.string_extract(filename, angle_min, start_string = angle_min_start, end_string = ',')
print angle_min

# Supreme angle for rotation if exists
angle_sup = [[] for n in filename]
angle_sup_start = "angle_sup="
angle_sup_label = "angle="

Label.string_extract(filename, angle_sup, start_string = angle_sup_start, end_string = ',')
print angle_sup

# Coupling strength if exists
coupling = [[] for n in filename]
coupling_start = "J="
coupling_label = "J="

Label.string_extract(filename, coupling, start_string = coupling_start, end_string = ',')
print coupling

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(label, length_label, length, end = True)
Label.label_build(label, coupling_label, coupling, end = True)
#Label.label_build(label, angle_min_label, angle_min, end = True)
Label.label_build(label, angle_sup_label, angle_sup, end = True)
#Label.label_build(label, realization_label, realization, end = True)
print label

matrix = []

index = 0
for n in filename:
    print "Read file", n
    matrix.append([])

    if cal_type[index] == "norm" and int(realization[index]) == 1:
        matrix[index] = read_matrix_norm(n)

    index += 1


matrix_off = []
for i in range(len(matrix)):
    matrix_off.append([])
    if len(matrix[i]) != 0:
        for row in range(len(matrix[i])):
            for col in range(len(matrix[i][row])):
                if row != col:
                    matrix_off[i].append(matrix[i][row][col])


import pylab

draw = []
index = 0
for i in range(len(matrix_off)):
    if len(matrix_off[i]) != 0:
        draw.append(Draw.Draw())

        draw[index].figure_init()
        draw[index].figure_set()

        bin_start = 0
        if coupling[i] != -1 and float(coupling[i])<0.5:
            bin_end = 0.01
            bin_width = 0.00005
        else:
            bin_end = 0.25
            bin_width = 0.005

        print label[i], min(matrix_off[i]), max(matrix_off[i])
        draw[index].hist(matrix_off[i], bin_start, bin_end, bin_width, label = label[i])
        pylab.legend(loc='upper right', ncol=1, prop={'size':15} )

        if coupling[i] != -1:
            J_s = coupling[i].replace('.','_')
        if angle_sup[i] != -1:
            angle_s = angle_min[i].replace('.','_')
        #pylab.xlabel(r"$\langle\Delta \phi_i \rangle$")

        pylab.subplots_adjust(bottom=0.12)

        print_label = label[i].replace('.', '_')
        print_label = print_label.replace(' ', '_')
        print_label = print_label + "_run=" + realization[i] + "_sigma_z_norm_bin"
        print print_label

        pylab.savefig(print_label+".pdf", box_inches='tight')

pylab.show()

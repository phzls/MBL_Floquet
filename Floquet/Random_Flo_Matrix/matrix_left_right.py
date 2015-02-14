path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

"""
This file plot relations between norms of matrix elements from leftmost as well as
rightmost sigma z matrices versus phases differences.
"""

import drawing as Draw
import label_build as Label
import numpy as np
from math import sqrt, pi

bin_num = 50
threshold = 0.03

def read_matrix_norm(filename):
    """ Read raw matrix elements norm"""
    matrix = np.loadtxt(filename+".txt")
    return matrix

def read_phases(filename):
    """ Read phases"""
    phase = np.loadtxt(filename+".txt")
    return phase

filename = [] # norm or phases can be added to it
label = []

right_norm_string = ",rightmost_sigma_z_norm_same"
left_norm_string = ",leftmost_sigma_z_norm_same"
phase_string = ",eval_phases_same"

# Obtain filenames
Label.name_read("left_right_name", filename, label)

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

# Realizations
realization = [[] for n in filename]
realization_string = "Realizations="
realization_label = "Run="

Label.string_extract(filename, realization, start_string = realization_string)

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
print label

Label.label_build(legend, '', general, end = True, concat = '')
Label.label_build(legend, length_label, length, end = True)
Label.label_build(legend, coupling_label, coupling, end = True)
#Label.label_build(label, angle_min_label, angle_min, end = True)
Label.label_build(legend, angle_sup_label, angle_sup, end = True)
Label.label_build(legend, realization_label, realization, end = True)

temp_left_matrix = []
temp_right_matrix = []
temp_prod_matrix = []
temp_phase = []

# Read files
for i in range(len(filename)):
    print "Read file", filename[i]

    temp_left_matrix.append([])
    temp_right_matrix.append([])
    temp_prod_matrix.append([])

    temp_left_matrix[i] = read_matrix_norm(filename[i]+left_norm_string)
    temp_right_matrix[i] = read_matrix_norm(filename[i]+right_norm_string)
    temp_prod_matrix[i] = np.multiply(temp_left_matrix[i], temp_right_matrix[i])

    temp_phase.append([])
    temp_phase[i] = read_phases(filename[i]+phase_string)

left_matrix = [] # Put matrix in an array. Maybe there is a native way to do this
right_matrix = []
prod_matrix = []
phase = []
phase_bin = [] # Binned phase
left_matrix_bin = [] # Average matrix elements in each bin
right_matrix_bin = []
prod_matrix_bin = []

bins = np.linspace(0, pi, bin_num+1)

# Construct data for plotting
for i in range(len(filename)):
    left_matrix.append([])
    right_matrix.append([])
    prod_matrix.append([])

    phase.append([])
    phase_bin.append([])
    left_matrix_bin.append([])
    right_matrix_bin.append([])
    prod_matrix_bin.append([])

    phase_bin[i] = [(bins[n]+bins[n+1])/2 for n in range(len(bins)-1)]
    left_matrix_bin[i] = [0 for n in phase_bin[i]]
    right_matrix_bin[i] = [0 for n in phase_bin[i]]
    prod_matrix_bin[i] = [0 for n in phase_bin[i]]
    count = [0 for n in phase_bin[i]] # Count number of matrix elements in each bin

    summation = 0

    print "Process", filename[i]

    for row in range(len(temp_left_matrix[i])):
        for col in range(len(temp_left_matrix[i][row])):
            left_matrix[i].append(temp_left_matrix[i][row][col])
            right_matrix[i].append(temp_right_matrix[i][row][col])
            prod_matrix[i].append(temp_prod_matrix[i][row][col])

            phase_diff = abs(temp_phase[i][row] - temp_phase[i][col])
            if phase_diff > pi:
                phase_diff = 2*pi - phase_diff
            phase[i].append(phase_diff)

            pos = 1
            while bins[pos] < phase_diff:
                pos += 1
            left_matrix_bin[i][pos-1] += temp_left_matrix[i][row][col]
            right_matrix_bin[i][pos-1] += temp_right_matrix[i][row][col]
            prod_matrix_bin[i][pos-1] += temp_prod_matrix[i][row][col]
            count[pos-1] += 1

    left_matrix_bin[i] = [left_matrix_bin[i][n]/float(count[n]) for n in
        range(len(left_matrix_bin[i]))]
    right_matrix_bin[i] = [right_matrix_bin[i][n]/float(count[n]) for n in
        range(len(left_matrix_bin[i]))]
    prod_matrix_bin[i] = [prod_matrix_bin[i][n]/float(count[n]) for n in
        range(len(left_matrix_bin[i]))]

    print sum(count)

import pylab
draw1 = Draw.Draw()

draw1.figure_init(width = 0)
draw1.figure_set()


must_in_range = ["J=0.1"]
must_not_range = []

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw1.plot(phase_bin, prod_matrix_bin, label = legend, trun = False)

pylab.legend(loc='center right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$\langle|\sigma^{z,l\times r}_{i,j}|\rangle$")
pylab.xlabel(r"$|\Delta\phi|$")

# If only one matrix is plotted
if len(draw1._plot_range) == 1:
    print_label = legend[draw1._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_bin_num_" + str(bin_num) + "_binned_ave_left_right_prod"

    print print_label

pylab.subplots_adjust(left=0.15)

pylab.savefig(print_label + ".pdf", box_inches='tight')

draw2 = Draw.Draw()

draw2.figure_init(width = 0)
draw2.figure_set()


must_in_range = ["J=0.1"]
must_not_range = []

draw2.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw2.plot(left_matrix_bin, right_matrix_bin, label = legend, trun = False)

pylab.legend(loc='center right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$\langle|\sigma^{z,r}_{i,j}|\rangle$")
pylab.xlabel(r"$\langle|\sigma^{z,l}_{i,j}|\rangle$")

# If only one matrix is plotted
if len(draw2._plot_range) == 1:
    print_label = legend[draw2._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_bin_num_" + str(bin_num) + "_binned_ave_left_vs_right"

    print print_label

pylab.subplots_adjust(left=0.15)
pylab.subplots_adjust(bottom=0.12)

pylab.savefig(print_label + ".pdf", box_inches='tight')

draw3 = Draw.Draw()

draw3.figure_init(width = 0)
draw3.figure_set()


must_in_range = ["J=0.1"]
must_not_range = []

draw3.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw3.plot(phase, prod_matrix, label = legend, trun = False)

pylab.legend(loc='center right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,l\times r}_{i,j}|$")
pylab.xlabel(r"$|\Delta\phi|$")

# If only one matrix is plotted
if len(draw2._plot_range) == 1:
    print_label = legend[draw2._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_left_right_prod"

    print print_label

pylab.savefig(print_label + ".png", box_inches='tight')

draw4 = Draw.Draw()

draw4.figure_init(width = 0)
draw4.figure_set()


must_in_range = ["J=0.1"]
must_not_range = []

draw4.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw4.plot(left_matrix, right_matrix, label = legend, trun = False)

pylab.legend(loc='center right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,r}_{i,j}|$")
pylab.xlabel(r"$|\sigma^{z,l}_{i,j}|$")

# If only one matrix is plotted
if len(draw2._plot_range) == 1:
    print_label = legend[draw2._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_left_vs_right"

    print print_label

pylab.savefig(print_label + ".png", box_inches='tight')

pylab.subplots_adjust(bottom=0.12)

draw5 = Draw.Draw()

draw5.figure_init(width = 0)
draw5.figure_set()


must_in_range = ["J=0.1"]
must_not_range = []

draw5.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw5.plot(phase_bin, right_matrix_bin, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,r}_{i,j}|$")
pylab.xlabel(r"$|\Delta\phi|$")

draw6 = Draw.Draw()

draw6.figure_init(width = 0)
draw6.figure_set()


must_in_range = ["J=0.1"]
must_not_range = []

draw6.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw6.plot(phase_bin, left_matrix_bin, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,l}_{i,j}|$")
pylab.xlabel(r"$|\Delta\phi|$")


pylab.show()

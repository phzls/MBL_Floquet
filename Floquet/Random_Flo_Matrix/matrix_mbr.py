path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

"""
This file plot norm of matrix elements versus phases differences.
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

norm_string = ",right_most_sigma_z_norm"
phase_string = ",eval_phases"

# Obtain filenames
Label.name_read("mbr_name", filename, label)

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

temp_matrix = []
temp_phase = []

# Read files
for i in range(len(filename)):
    print "Read file", filename[i]

    temp_matrix.append([])
    temp_matrix[i] = read_matrix_norm(filename[i]+norm_string)

    temp_phase.append([])
    temp_phase[i] = read_phases(filename[i]+phase_string)

matrix = []
phase = []
phase_bin = [] # Binned phase
matrix_bin = [] # Average matrix elements in each bin
matrix_count = [] # Number of matrix elments above a threshold in each bin

bins = np.linspace(0, pi, bin_num+1)

# Construct data for plotting
for i in range(len(filename)):
    matrix.append([])
    phase.append([])
    phase_bin.append([])
    matrix_bin.append([])
    matrix_count.append([])

    phase_bin[i] = [(bins[n]+bins[n+1])/2 for n in range(len(bins)-1)]
    matrix_count[i] = [0 for n in phase_bin[i]]
    matrix_bin[i] = [0 for n in phase_bin[i]]
    count = [0 for n in phase_bin[i]] # Count number of matrix elements in each bin

    summation = 0

    print "Process", filename[i]

    for row in range(len(temp_matrix[i])):
        for col in range(len(temp_matrix[i][row])):
            matrix[i].append(temp_matrix[i][row][col])

            summation += temp_matrix[i][row][col] ** 2

            phase_diff = abs(temp_phase[i][row] - temp_phase[i][col])
            if phase_diff > pi:
                phase_diff = 2*pi - phase_diff
            phase[i].append(phase_diff)

            pos = 1
            while bins[pos] < phase_diff:
                pos += 1
            matrix_bin[i][pos-1] += temp_matrix[i][row][col]
            count[pos-1] += 1

            if temp_matrix[i][row][col] > threshold:
                matrix_count[i][pos-1] += 1

    matrix_bin[i] = [matrix_bin[i][n]/float(count[n]) for n in range(len(matrix_bin[i]))]


    print filename[i], summation
    print len(phase_bin[i]), len(matrix_bin[i]), len(matrix_count[i])

import pylab
draw1 = Draw.Draw()

draw1.figure_init(width = 0)
draw1.figure_set()


must_in_range = ["J=0.4", "Run=1", "Inter"]
must_not_range = []

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw1.plot(phase, matrix, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^z_{i,j}|$")
pylab.xlabel(r"$|\Delta\phi|$")

# If only one matrix is plotted
if len(draw1._plot_range) == 1:
    print_label = legend[draw1._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_MBR"

#pylab.savefig(print_label + ".png", box_inches='tight')

draw2 = Draw.Draw()

draw2.figure_init(width = 0)
draw2.figure_set()


#must_in_range = ["J=0.1", "Run=1"]
#must_not_range = ["Inter"]

draw2.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw2.plot(phase_bin, matrix_bin, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$\langle|\sigma^z_{i,j}|\rangle$")
pylab.xlabel(r"$|\Delta\phi|$")

# If only one matrix is plotted
if len(draw2._plot_range) == 1:
    print_label = legend[draw2._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_bin_num_" + str(bin_num) + "_binned_ave_MBR"

#pylab.savefig(print_label + ".pdf", box_inches='tight')

draw3 = Draw.Draw()

draw3.figure_init(width = 0)
draw3.figure_set()


#must_in_range = ["J=0.1", "Run=1"]
#must_not_range = ["Inter"]

draw3.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw3.plot(phase_bin, matrix_count, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel("Count")
pylab.xlabel(r"$|\Delta\phi|$")

threshold_plot = str(threshold).replace('.','_')

# If only one matrix is plotted
if len(draw3._plot_range) == 1:
    print_label = legend[draw3._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_bin_num_" + str(bin_num) + "_threshold_" + threshold_plot + "_count_MBR"

#pylab.savefig(print_label + ".pdf", box_inches='tight')

pylab.show()

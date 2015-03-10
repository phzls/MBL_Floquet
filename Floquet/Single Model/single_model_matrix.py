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
import networkx as nx
from math import sqrt, pi, floor, log

# single_model file name index
name_index = 2

bin_num = 20

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

right_norm_string = ",rightmost_sigma_z_under_one_model_norm,v2"
left_norm_string = ",leftmost_sigma_z_under_one_model_norm,v2"
phase_string = ",eval_phases_under_one_model,v2"

# Obtain filenames
Label.name_read("single_model_name"+str(name_index), filename, label, exclude_s = "Markov")

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
#Label.label_build(legend, coupling_label, coupling, end = True)
#Label.label_build(label, angle_min_label, angle_min, end = True)
Label.label_build(legend, angle_sup_label, angle_sup, end = True)
#Label.label_build(legend, realization_label, realization, end = True)

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
bin_width = bins[1]-bins[0]

# Construct data for plotting
# Diagonal terms are excluded
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
            if row == col:
                continue

            left_matrix[i].append(temp_left_matrix[i][row][col])
            right_matrix[i].append(temp_right_matrix[i][row][col])
            prod_matrix[i].append(temp_prod_matrix[i][row][col])

            phase_diff = abs(temp_phase[i][row] - temp_phase[i][col])
            if phase_diff > pi:
                phase_diff = 2*pi - phase_diff
            phase[i].append(phase_diff)

            pos = int(floor(phase_diff/bin_width))

            left_matrix_bin[i][pos] += temp_left_matrix[i][row][col]
            right_matrix_bin[i][pos] += temp_right_matrix[i][row][col]
            prod_matrix_bin[i][pos] += temp_prod_matrix[i][row][col]
            count[pos] += 1

    left_matrix_bin[i] = [left_matrix_bin[i][n]/float(count[n]) for n in
        range(len(left_matrix_bin[i]))]
    right_matrix_bin[i] = [right_matrix_bin[i][n]/float(count[n]) for n in
        range(len(left_matrix_bin[i]))]
    prod_matrix_bin[i] = [prod_matrix_bin[i][n]/float(count[n]) for n in
        range(len(left_matrix_bin[i]))]

    print sum(count)

    print "In the first bin: ", count[0]
    print "In the middle of prod matrix bin: ", prod_matrix_bin[i][10]
    print "Smallest in the prod matrix bin", min(prod_matrix_bin[i])

# Construct an adjacency matrix with uniform weight 1 for largest off-diagonal 
# element in each row

adjacency = []
max_value = []
max_value_phase_diff = []
in_degrees = []
for i in range(len(filename)):
    row = len(temp_prod_matrix[i])
    adjacency.append(np.zeros((row,row)))
    max_value.append([])
    max_value_phase_diff.append([])
    in_degrees.append(np.zeros(row))

    prod = 0.0

    np.fill_diagonal(temp_prod_matrix[i],0)
    index_array = np.argmax(temp_prod_matrix[i], axis = 1)

    for n in range(row):
        m = index_array[n]
        adjacency[-1][n][m] = 1

        max_value[-1].append(temp_prod_matrix[i][n][m])

        delta_phi = abs(temp_phase[i][m]-temp_phase[i][n])
        if delta_phi > pi:
            delta_phi = 2*pi - delta_phi
        max_value_phase_diff[-1].append(delta_phi)

        in_degrees[-1][m] += 1

        prod += log(temp_prod_matrix[i][n][m])
    print "Log Product of maximum values:", prod
    print "Maximum: ", max(max_value[-1])

# Construct graphs
G = []
for i in range(len(filename)):
    G.append( nx.from_numpy_matrix(adjacency[i],nx.DiGraph()) )


import pylab

"""
draw1 = Draw.Draw()

draw1.figure_init(width = 0, logy = True)
draw1.figure_set()


must_in_range = ["J=0.3"]
must_not_range = []

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw1.plot(phase_bin, prod_matrix_bin, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$\langle|\sigma^{z,l\times r}_{i,j}|\rangle$")
pylab.xlabel(r"$|\Delta\phi|$")

# If only one matrix is plotted
if len(draw1._plot_range) == 1:
    print_label = legend[draw1._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_bin_num_" + str(bin_num) + "_binned_ave_left_right_prod_log"

    print print_label

pylab.subplots_adjust(left=0.15)

pylab.savefig(print_label + ",v2.pdf", box_inches='tight')

draw2 = Draw.Draw()

draw2.figure_init(width = 0,logy = True, logx = True)
draw2.figure_set()


must_in_range = ["J=0.3"]
must_not_range = []

draw2.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw2.plot(left_matrix_bin, right_matrix_bin, label = legend, trun = False)

pylab.legend(loc='upper left', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$\langle|\sigma^{z,r}_{i,j}|\rangle$")
pylab.xlabel(r"$\langle|\sigma^{z,l}_{i,j}|\rangle$")

# If only one matrix is plotted
if len(draw2._plot_range) == 1:
    print_label = legend[draw2._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_bin_num_" + str(bin_num) + "_binned_ave_left_vs_right_loglog"

    print print_label

pylab.subplots_adjust(left=0.13)
pylab.subplots_adjust(bottom=0.13)

pylab.savefig(print_label + ",v2.pdf", box_inches='tight')

draw3 = Draw.Draw()

draw3.figure_init(width = 0,logy = True)
draw3.figure_set()


must_in_range = ["J=0.3"]
must_not_range = []

draw3.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw3.plot(phase, prod_matrix, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,l\times r}_{i,j}|$")
pylab.xlabel(r"$|\Delta\phi|$")

# If only one matrix is plotted
if len(draw2._plot_range) == 1:
    print_label = legend[draw2._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_left_right_prod_log"

    print print_label

pylab.subplots_adjust(left=0.14)
pylab.savefig(print_label + ",v2.png", box_inches='tight')

draw4 = Draw.Draw()

draw4.figure_init(width = 0,logy = True, logx = True)
draw4.figure_set()


must_in_range = ["J=0.3"]
must_not_range = []

draw4.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw4.plot(left_matrix, right_matrix, label = legend, trun = False)

pylab.legend(loc='lower right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,r}_{i,j}|$")
pylab.xlabel(r"$|\sigma^{z,l}_{i,j}|$")

# If only one matrix is plotted
if len(draw2._plot_range) == 1:
    print_label = legend[draw2._plot_range[0]].replace('.', '_')
    print_label = print_label.replace(' ', '_')
    print_label = print_label + "_left_vs_right_log"

    print print_label

pylab.subplots_adjust(left=0.14)
pylab.subplots_adjust(bottom=0.14)
pylab.savefig(print_label + ",v2.png", box_inches='tight')

pylab.subplots_adjust(bottom=0.12)

draw5 = Draw.Draw()

draw5.figure_init(width = 0,logy = True)
draw5.figure_set()


must_in_range = ["J=0.3"]
must_not_range = []

draw5.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw5.plot(phase_bin, right_matrix_bin, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,r}_{i,j}|$")
pylab.xlabel(r"$|\Delta\phi|$")
pylab.subplots_adjust(left=0.15)

draw6 = Draw.Draw()

draw6.figure_init(width = 0,)
draw6.figure_set()


must_in_range = ["J=0.3"]
must_not_range = []

draw6.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw6.plot(phase_bin, left_matrix_bin, label = legend, trun = False)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))
pylab.ylabel(r"$|\sigma^{z,l}_{i,j}|$")
pylab.xlabel(r"$|\Delta\phi|$")
pylab.subplots_adjust(left=0.15)

# Plot histogram and visualize graphs for one set of data
data_index = 0

draw7 = Draw.Draw()

draw7.figure_init()
draw7.figure_set()

draw7.hist(max_value[data_index], bin_num = 30, label = "maximum values")

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))


print_label = legend[data_index].replace('.', '_')
print_label = print_label.replace(' ', '_')
print_label = print_label + "_maximum_value_per_row_hist"

print print_label

pylab.savefig(print_label + ",v2.pdf", box_inches='tight')

draw8 = Draw.Draw()

draw8.figure_init()
draw8.figure_set()

draw8.hist(max_value_phase_diff[data_index], bin_num = 30, label = "phases")

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))

print_label = legend[data_index].replace('.', '_')
print_label = print_label.replace(' ', '_')
print_label = print_label + "_maximum_value_per_row_phase_diff_abs_hist"

print print_label

pylab.savefig(print_label + ",v2.pdf", box_inches='tight')

pylab.figure(9)
nx.draw(G[data_index])

print_label = legend[data_index].replace('.', '_')
print_label = print_label.replace(' ', '_')
print_label = print_label + "_maximum_value_per_row_uniform_adjacency_graph"

print print_label

pylab.savefig(print_label + ",v2.pdf", box_inches='tight')



draw10 = Draw.Draw()

draw10.figure_init()
draw10.figure_set(fig_num = 10)

hist, bin_edges = np.histogram(in_degrees[0])

print bin_edges
print hist

draw10.hist(in_degrees[0], label = "In Degrees")
pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.1, 1))

print_label = legend[0].replace('.', '_')
print_label = print_label.replace(' ', '_')
print_label = print_label + "_maximum_value_per_row_in_degree_hist"

print print_label

pylab.savefig(print_label + ",v1.pdf", box_inches='tight')

pylab.show()

"""
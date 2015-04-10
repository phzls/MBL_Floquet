__author__ = 'liangshengzhang'

path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing_new_1 as Draw
import label_build as Label
import numpy as np

include_s = "component" # A string that must be included in filename construct
exclude_s = "Shift" # A string that must be excluded in filename construct
ver_num = 1
scale = False
scale_target = 0.48

def read_file(filename, realization = 1):
    """  Read averages and standard deviation (if exist) from files
         concerning transition of floquet operators."""

    temp = np.loadtxt(filename+".txt", skiprows=1)

    num_parameter = len(temp[0]) - 1
    if num_parameter % 2 != 0:
        raise Exception("Not all parameters have both mean and sd.")

    W = []
    ave = [[] for n in range(num_parameter/2)]
    err = [[] for n in range(num_parameter/2)]

    for n in range(len(temp)):
        W.append(temp[n][0])
        for i in range(num_parameter/2):
            ave[i].append(temp[n][2*i+1])
            if len(temp[n]) > 2:
                err[i].append(temp[n][2*i+2] / (realization**0.5) )
    return W, ave, err

def scale_point(data, target):
    """
    This function finds the point where the data first becomes no less than the target.
    The data is assumed to be in ascending order. If not found, none will be returned.
    :param data: the data file
    :param target: the target value to be not less than
    :return: the position in the data. If not found, none will be returned
    """
    L =len(data)
    pos = 0
    while data[pos] < target:
        pos += 1
        if pos >= L:
            pos = None
            break
    return pos

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

# Realizations
realization = [[] for n in filename]
realization_string = "Run="
realization_label = "Run="

Label.string_extract(filename, realization, start_string = realization_string, end_string = ',')

# Number of points
J_N = [[] for n in filename]
J_N_string = "J_N="

Label.string_extract(filename, J_N, start_string = J_N_string, end_string = ',')

# Version number if exists
version = [[] for n in filename]
version_start = ",v"
version_label = "v"

Label.string_extract(filename, version, start_string = version_start)

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(legend, '', general, end = True, concat = '')

Label.label_build(label, length_label, length, end = True)
Label.label_build(legend, length_label, length, end = True)

Label.label_build(label, realization_label, realization, end = True)
Label.label_build(label, version_label, version, end = True)

#print label
#print legend

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

W = []
ave = []
err = []

index = 0
for n in filename:
    print "Read file", n
    temp_W, temp_ave, temp_err = read_file(filename[index], float(realization[index]))
    W.append(temp_W)
    ave.append(temp_ave)
    err.append(temp_err)
    index += 1

"""
scaled_data = [[],[],[]]
for n in range(len(data[0])):
    scale_pos = scale_point(data[0][n], scale_target)
    scale_factor = data[1][n][scale_pos]
    scaled_data[0].append([x for x in data[0][n]])
    scaled_data[1].append([x for x in data[1][n]/scale_factor])
    scaled_data[2].append([x for x in data[2][n]/scale_factor])
"""

import pylab

L = 12

log_y = False

# Find the index corresponding to this L
for i in range(len(filename)):
    if int(length[i]) == L:
        L_instance = i
        break

draw1 = Draw.Draw()

draw1.plot_set(logy=log_y)

must_in_range = ["v"+str(ver_num)]
must_not_range = []

plot_range = draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range)

legend = ["<zz>^2", "<z>^2<z>^2", "(<z><z>)^2"]

W_temp = []
for i in range(len(ave[L_instance])):
    W_temp.append(W[L_instance])

plot_range = range(len(ave[L_instance]))

draw1.errorbar(W_temp, ave[L_instance], yerr = err[L_instance], plot_range= plot_range,
               legend = legend)

pylab.legend(loc='upper left', ncol=1, prop={'size':14})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel("Time Correlation Components")
pylab.xlabel("W")

num_pts = J_N[draw1._plot_range[0]]
num_run = realization[draw1._plot_range[0]]

save_name = "Random_Simple_Floquet_" + "L_" + str(L) + "_J_N_" + str(num_pts) + "_Run_" \
            + str(num_run) + "_" + "zz_time_corr_component"

"""
if scale:
    save_name += "_scaled_target_" + str(scale_target).replace('.','_')
"""
if draw1._xmax < max(W[L_instance]):
    save_name += "_zoom"


if log_y:
    save_name += "_logy"
print save_name

#pylab.subplots_adjust(left=0.14)
pylab.savefig(save_name + "_v" + str(ver_num) + ".pdf", box_inches='tight')


pylab.show()


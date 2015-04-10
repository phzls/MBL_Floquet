__author__ = 'liangshengzhang'

"""
Properties at single J
"""

path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing_new_1 as Draw
import label_build as Label
import numpy as np

include_s = None # A string that must be included in filename construct
exclude_s = None # A string that must be excluded in filename construct
ver_num = 1

def read_file(filename, L, distance, ave, err, realization = 1):
    """  Read averages and standard deviation (if exist) from files
         concerning transition of floquet operators."""

    temp = np.loadtxt(filename+".txt")

    for n in range(len(temp)):
        distance.append(L - 2*(temp[n][0]+1))
        ave.append(temp[n][1])
        if len(temp[n]) > 2:
            err.append(temp[n][2] / (realization**0.5) )

filename = []
label = [] # General label for excuding when plotting

# Obtain filenames
Label.name_read("single_name", filename, label, include_s = include_s, exclude_s=exclude_s)

legend = ['' for n in filename] # labels used as legends in plotting

# Generic labels
general = [[] for n in filename]

Label.string_extract(filename, general, end_string = ",size=")
Label.label_change(general)

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

# J value
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
Label.label_build(label, realization_label, realization, end = True)
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

data = [[],[],[]]

index = 0
for n in filename:
    print "Read file", n
    for i in range(3):
        data[i].append([])
    L = float(length[index])
    print "Length: ", L
    read_file(filename[index], L, data[0][-1], data[1][-1], data[2][-1], float(realization[index]))
    index += 1

import pylab
draw1 = Draw.Draw()

draw1.plot_set()

must_in_range = ["v"+str(ver_num)]
must_not_range = []

plot_range = draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range)
#print plot_range

plot_data = data

if len(data[2][0]) > 0:
    draw1.errorbar(plot_data[0], plot_data[1], yerr = plot_data[2], legend = legend)
else:
    draw1.plot(plot_data[0], plot_data[1], legend = legend)

pylab.legend(loc='upper right', ncol=1, prop={'size':14})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel("ZZ Correlation")
pylab.xlabel("Distance")

J_val = str(J[plot_range[0]]).replace('.','_')
num_run = realization[draw1._plot_range[0]]

#save_name = "Random_Simple_Floquet_" + "J_N_" + str(num_pts) + "_Run_" + str(num_run) + "_" + "zz_corr"
save_name = "Random_Simple_Floquet_zz_all_corr_" + "J_" + J_val

print save_name

#pylab.subplots_adjust(left=0.14)
pylab.savefig(save_name + "_v" + str(ver_num) + ".pdf",box_inches='tight')

pylab.show()
__author__ = 'liangshengzhang'

path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing as Draw
import label_build as Label
import numpy as np

include_s = None # A string that must be included in filename construct
excludes_s = None # A string that must be excluded in filename construct
ver_num = 1
scale = False

def read_file(filename, W, ave, err, realization = 1):
    """  Read averages and standard deviation (if exist) from files
         concerning transition of floquet operators."""

    temp = np.loadtxt(filename+".txt")

    for n in range(len(temp)):
        W.append(temp[n][0])
        ave.append(temp[n][1])
        if len(temp[n]) > 2:
            err.append(temp[n][2] / (realization**0.5) )

filename = []
label = [] # General label for excuding when plotting

# Obtain filenames
Label.name_read("name", filename, label, include_s = include_s, exclude_s = excludes_s)

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

# Type considered
type = [[] for n in filename]
type_start = "zz_"
type_label = "Type="

Label.string_extract(filename, type, start_string= type_start)

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(legend, '', general, end = True, concat = '')

Label.label_build(label, length_label, length, end = True)
Label.label_build(legend, length_label, length, end = True)

Label.label_build(label, realization_label, realization, end = True)
Label.label_build(label, version_label, version, end = True)
Label.label_build(label, type_label, type, end = True)

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

#print label
#print legend

data = [[],[],[]]

index = 0
for n in filename:
    print "Read file", n
    for i in range(3):
        data[i].append([])
    read_file(filename[index], data[0][-1], data[1][-1], data[2][-1], float(realization[index]))
    index += 1

scaled_data = [[],[],[]]
for n in range(len(data[0])):
    scale_factor = data[1][n][0] # Use average at W=0 as scaling factor
    scaled_data[0].append([x for x in data[0][n]])
    scaled_data[1].append([x for x in data[1][n]/scale_factor])
    scaled_data[2].append([x for x in data[2][n]/scale_factor])

import pylab
draw1 = Draw.Draw()

draw1.figure_init()
draw1.figure_set()

must_in_range = ["v"+str(ver_num), "L=10", "Type"]
must_not_range = ["time"]

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

plot_data = data
if scale:
    plot_data = scaled_data

if len(data[2][0]) > 0:
    draw1.errorbar(plot_data[0], plot_data[1], yerr = plot_data[2], label = legend)
else:
    draw1.plot(plot_data[0], plot_data[1], label=legend)

pylab.legend(loc='lower left', ncol=1, prop={'size':14})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel("ZZ Correlation")
pylab.xlabel("W")

num_pts = J_N[draw1._plot_range[0]]
num_run = realization[draw1._plot_range[0]]

save_name = "Random_Floquet_" + "J_N_" + str(num_pts) + "_Run_" + str(num_run) + "_" + "zz_corr_compare"

if scale:
    save_name += "_scaled"

print save_name

pylab.subplots_adjust(left=0.14)
pylab.savefig(save_name + "_v" + str(ver_num) + ".pdf",box_inches='tight')

pylab.show()


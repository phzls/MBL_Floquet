path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing as Draw
import label_build as Label
import numpy as np

def read_mean(filename, J, ave, err, realization = 1):
    """  Read averages and standard deviation of averages from files
         concerning mean values of level statistics."""

    temp = np.loadtxt(filename+".txt")

    for n in range(len(temp)):
        J.append(temp[n][0])
        ave.append(temp[n][1])
        err.append(temp[n][2] / (realization**0.5) )



#    f = open(filename+".txt")
#    for line in f:
#        J.append(float(line.split()[0]))
#        ave.append(float(line.split()[1]))
#        err.append(float(line.split()[2]) / (realization**0.5) )


def read_level(filename, J, level):
    """ Read raw level statistitcs"""

    temp = np.loadtxt(filename+".txt")

    for n in range(len(temp)):
        J.append(temp[n][0])
        level.append([])
        level[-1] = temp[n][1:]

#    f = open(filename+".txt")
#    for line in f:
#        J.append(float(line.split()[0]))
#        level.append([])
#        for i in xrange(len(line.split())-1):
#            level[-1].append(float(line.split()[i+1]))



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
type_string = "level"

Label.string_extract(filename, cal_type, desired_string = type_string)
Label.label_change(cal_type)

# Realizations
realization = [[] for n in filename]
realization_string = "Realizations="

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
angle_sup_label = "angle sup="

Label.string_extract(filename, angle_sup, start_string = angle_sup_start, end_string = ',')
print angle_sup

# Number of points
J_N = [[] for n in filename]
J_N_string = "J_N="

Label.string_extract(filename, J_N, start_string = J_N_string, end_string = ',')

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(label, '', cal_type, end = True)
Label.label_build(label, length_label, length, end = True)
Label.label_build(label, angle_min_label, angle_min, end = True)
Label.label_build(label, angle_sup_label, angle_sup, end = True)
Label.label_build(label, "runs=", realization, end = True)
Label.label_build(label, J_N_string, J_N, end = True)
print label


data = [[],[],[]]

index = 0
for n in filename:
    print "Read file", n
    for i in range(3):
        data[i].append([])
    if n.find("mean") > -1:
        read_mean(filename[index], data[0][-1], data[1][-1], data[2][-1],
                 float(realization[index]))
    else:
        read_level(filename[index], data[0][-1], data[1][-1])

    index += 1



import pylab
draw1 = Draw.Draw()

draw1.figure_init(ymax = 2.1, ymin = 1.2)
draw1.figure_set()

in_range = ["L=11"] # Angles or length
must_in_range = ["square", "runs=100"]

extra_index= [] # The original two scenarios
"""
for n in xrange(len(label)):
    if label[n].find("square") > -1 and label[n].find("angle") == -1:
        extra_index.append(n)
"""

plot_label = ['' for n in filename]
print general
Label.label_build(plot_label, '', general, end = True, concat = '')
Label.label_build(plot_label, length_label, length, end = True)
Label.label_build(plot_label, "angle=", angle_min, end = True)
Label.label_build(plot_label, "runs=", realization, end = True)

print "Plot label:", plot_label

draw1.plot_range(label, in_range = in_range, must_in_range = must_in_range,
                 index = extra_index, printout = True)

draw1.errorbar(data[0], data[1], yerr = data[2], label = plot_label)

pylab.legend(loc='upper right', ncol=1, prop={'size':15})#, bbox_to_anchor=(1.15, 0.8))
pylab.ylabel(r"$\langle(\Delta\phi)^2\rangle$")
pylab.xlabel(r"$j$")

#pylab.savefig("Inter_Floquet_12_level_spacing_mean_square_compare_v1.pdf",box_inches='tight')


draw2 = Draw.Draw()

draw2.figure_init()
draw2.figure_set()

index = 0
for n in filename:
    if n.find("mean") == -1 and n.find("Inter") > -1:
        break
    else:
        index += 1

print index, len(data[0])

print min(data[1][index][1]), max(data[1][index][1])
print min(data[1][index][9]), max(data[1][index][9])

bin_width = 0.05

instance = 1

label = ( general[index] + " L=" + length[index]+" "+"J="+str(data[0][index][instance])
          )
print label

draw2.hist(data[1][index][instance], bin_width = bin_width, label = label)
pylab.legend(loc='upper right', ncol=1, prop={'size':15} )

J_s = str(data[0][index][instance]).replace('.','_')
#angle_s = angle_min[index].replace('.','_')
pylab.xlabel(r"$\langle\Delta \phi_i \rangle$")

print_label = label.replace('.', '_')
print_label = print_label.replace(' ', '_')
print_label = print_label + "_level_spacing"

print print_label

pylab.subplots_adjust(bottom=0.12)

#pylab.savefig(print_label+".pdf", box_inches='tight')

pylab.show()





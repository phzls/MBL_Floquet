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

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(label, '', cal_type, end = True)
Label.label_build(label, length_label, length, end = True)
print label

# Realizations
realization = [[] for n in filename]
realization_string = "Realizations="

Label.string_extract(filename, realization, start_string = realization_string, end_string = ',')
print realization

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

draw1.figure_init()
draw1.figure_set()

in_range = ["square"]

plot_label = ['' for n in filename]
print general
Label.label_build(plot_label, '', general, end = True, concat = '')
Label.label_build(plot_label, length_label, length, end = True)
print plot_label

draw1.plot_range(label, in_range = in_range)

draw1.errorbar(data[0], data[1], yerr = data[2], label = plot_label)

pylab.legend(loc='upper right', ncol=1, prop={'size':15} )
pylab.ylabel(r"$\langle(\Delta\phi)^2\rangle$")
pylab.xlabel(r"$j$")

#pylab.savefig("Floquet_10_level_spacing_mean_square_compare.pdf",box_inches='tight')


draw2 = Draw.Draw()

draw2.figure_init()
draw2.figure_set()

index = 0
for n in filename:
    if n.find("mean") == -1 and n.find("Rotation") > -1:
        break
    else:
        index += 1

print index, len(data[0])

print min(data[1][index][1]), max(data[1][index][1])
print min(data[1][index][9]), max(data[1][index][9])

bin_width = 0.05

instance = 1

label = general[index] + " L=" + length[index]+" "+"J="+str(data[0][index][instance])
print label

draw2.hist(data[1][index][instance], bin_width, label = label)
pylab.legend(loc='upper right', ncol=1, prop={'size':15} )

J_s = str(data[0][index][instance]).replace('.','_')
pylab.xlabel(r"$\langle\Delta \phi_i \rangle$")

pylab.subplots_adjust(bottom=0.12)

pylab.savefig("Random_Rotation_Floquet_10_level_spacing_J_"+J_s+"_v2.pdf",
              box_inches='tight')

pylab.show()





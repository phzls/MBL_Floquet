path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing as Draw
import label_build as Label
import numpy as np

def read_ent(filename, t, ent_ave, ent_err, realization = 1):
    """  Read averages and standard deviation of entropies from files
         concerning time evolution of entropy per mode."""

    temp = np.loadtxt(filename+".txt")

    for n in range(len(temp)):
        t.append(temp[n][0])
        ent_ave.append(temp[n][1])
        ent_err.append(temp[n][2] / (realization**0.5) )

filename = []
label = [] # General label for excuding when plotting

# Obtain filenames
Label.name_read("name", filename, label)

legend = ['' for n in filename] # labels used as legends in plotting

# Generic labels
general = [[] for n in filename]

Label.string_extract(filename, general, end_string = "_L=")
Label.label_change(general)

#print general

# Length
length = [[] for n in filename]
length_start = "L="
length_end = ","
length_label = "L="

Label.string_extract(filename, length, start_string = length_start, end_string = length_end)

# Type of calculation
""" I should put a more easily recordnized label in front of it """
cal_type = [[] for n in filename]
type_string = "entropy"
type_label = "type is "

Label.string_extract(filename, cal_type, desired_string = type_string)
Label.label_change(cal_type)

# Realizations
realization = [[] for n in filename]
realization_string = "Realizations="
realization_label = "Run="

Label.string_extract(filename, realization, start_string = realization_string, end_string = ',')

# Initial conditions
init = [[] for n in filename]
init_string = "Init_"

Label.string_extract(filename, init, start_string = init_string, end_string = ',')
Label.label_change(init)

# Total time step
Total = [[] for n in filename]
Total_string = "Total_time_step="
Total_label = "Total="

Label.string_extract(filename, Total, start_string = Total_string, end_string = ',')

# jump
jump = [[] for n in filename]
jump_string = "jump="
jump_label = "jump="

Label.string_extract(filename, jump, start_string = jump_string, end_string = ',')

# minimum angle for rotation if exists
angle_min = [[] for n in filename]
angle_min_start = "angle_min="
angle_min_label = "angle min="

Label.string_extract(filename, angle_min, start_string = angle_min_start, end_string = ',')

# Supreme angle for rotation if exists
angle_sup = [[] for n in filename]
angle_sup_start = "angle_sup="
angle_sup_label = "angle sup="

Label.string_extract(filename, angle_sup, start_string = angle_sup_start, end_string = ',')

# Coupling strength if exists
coupling = [[] for n in filename]
coupling_start = "J="
coupling_label = "J="

Label.string_extract(filename, coupling, start_string = coupling_start, end_string = ',')

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(legend, '', general, end = True, concat = '')

Label.label_build(label, length_label, length, end = True)
Label.label_build(legend, length_label, length, end = True)

Label.label_build(label, coupling_label, coupling, end = True)
Label.label_build(legend, coupling_label, coupling, end = True)

Label.label_build(label, angle_min_label, angle_min, end = True)

for i in range(len(filename)):
    # If angle_min is not the same as angle_sup, assume it's full angle
    if angle_min[i] != -1:
        if angle_min[i] == angle_sup[i]:
            legend[i] += " angle=" + angle_min[i]
        else:
            legend[i] += " full angle"

Label.label_build(label, angle_sup_label, angle_sup, end = True)

Label.label_build(label, realization_label, realization, end = True)
Label.label_build(label, type_label, cal_type, end = True)
Label.label_build(label, Total_label, Total, end = True)
Label.label_build(label, jump_label, jump, end = True)

Label.label_build(label, '', init, end = True)
Label.label_build(legend, '', init, end = True)



#print label
#print legend

time = []
ent_ave = []
ent_err = []

for n in range(len(filename)):
    time.append([])
    ent_ave.append([])
    ent_err.append([])
    read_ent(filename[n], time[-1], ent_ave[-1], ent_err[-1],float(realization[n]))

print len(time[12]), len(ent_ave[12])

import pylab
draw1 = Draw.Draw()

draw1.figure_init(xmin=0.8, ymin = 0, ymax = 4, logx = True)
draw1.figure_set()

#in_range = ["3.14", "1.04", "1.57", "2.09", "2.61"] # Angles
must_in_range = ["Inter Random Floquet", "entropy", "Run=400 ", "Total=30 "]
must_not_range = ["J=0.9"]

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw1.errorbar(time, ent_ave, yerr = ent_err, label = legend)

pylab.legend(loc='upper left', ncol=1, prop={'size':14})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel(r"$S(t)$")
pylab.xlabel("time")

pylab.savefig("Inter_Random_Floquet_10_entropy_compare_long_log_time,v1.pdf",box_inches='tight')

pylab.show()
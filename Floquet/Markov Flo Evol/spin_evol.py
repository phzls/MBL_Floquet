path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing as Draw
import label_build as Label
import numpy as np

def read_spin(filename, t, spin_ave, spin_err, realization = 1, model = 1):
    """  Read averages and standard deviation of average spin values from files
         concerning time evolution of entropy per mode. If realization number is
         1, model number is used to calculate the sd of mean."""

    temp = np.loadtxt(filename+".txt")

    num = realization
    if num == 1 or num == -1:
        num = model

    for n in range(len(temp)):
        t.append(temp[n][0])
        spin_ave.append(temp[n][1])
        spin_err.append(temp[n][2] / (num**0.5) )

filename = []
label = [] # General label for excuding when plotting

# Obtain filenames
Label.name_read("spin_name", filename, label)

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
type_string = "leftmost"
type_label = "type is "

Label.string_extract(filename, cal_type, desired_string = type_string)
Label.label_change(cal_type)

# Task of calculation
cal_task = [[] for n in filename]
task_string = "Task_"
task_label = "task is "

Label.string_extract(filename, cal_task, start_string = task_string, end_string = ',')
Label.label_change(cal_task)

# Realizations
realization = [[] for n in filename]
realization_string = "Realizations="
realization_label = "Run="

Label.string_extract(filename, realization, start_string = realization_string, end_string = ',')

# Models
model_num = [[] for n in filename]
model_num_string = "Model_Num="
model_num_label = "Model_Num="

Label.string_extract(filename, model_num, start_string = model_num_string, end_string = ',')

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
#Label.label_build(legend, '', general, end = True, concat = '')

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

Label.label_build(label, model_num_label, model_num, end = True)
Label.label_build(legend, model_num_label, model_num, end = True)

Label.label_build(label, type_label, cal_type, end = True)
Label.label_build(label, task_label, cal_task, end = True)
Label.label_build(label, Total_label, Total, end = True)
Label.label_build(label, jump_label, jump, end = True)

Label.label_build(label, '', init, end = True)
Label.label_build(legend, '', init, end = True)



print label
#print legend

time = []
spin_ave = []
spin_err = []

for n in range(len(filename)):
    time.append([])
    spin_ave.append([])
    spin_err.append([])
    read_spin(filename[n], time[-1], spin_ave[-1], spin_err[-1],float(realization[n]),
        float(model_num[n]))


import pylab
draw1 = Draw.Draw()

draw1.figure_init()
draw1.figure_set()

#in_range = ["3.14", "1.04", "1.57", "2.09", "2.61"] # Angles
must_in_range = ["Markov Inter Random Both X Floquet", "L=8"]
must_not_range = ["J=0.1"]

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw1.errorbar(time, spin_ave, yerr = spin_err, label = legend)

pylab.legend(loc='left', ncol=1, prop={'size':14})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel(r"$\sigma^{z,l}(t)$")
pylab.xlabel("time")

#pylab.savefig("Inter_Random_Floquet_Markov_8_J_0_9_largest_left_spin_eigenstate_left_spin_compare_log.png",box_inches='tight')

pylab.show()

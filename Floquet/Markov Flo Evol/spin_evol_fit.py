path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing as Draw
import label_build as Label
import numpy as np
from scipy import stats
from math import sqrt, log

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
Label.name_read("spin_disorder_name", filename, label)
#Label.name_read("spin_thermal_name", filename, label)

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
jump_string = ",jump="
jump_label = "jump="

Label.string_extract(filename, jump, start_string = jump_string, end_string = ',')

# Markov time jump
markov_time_jump = [[] for n in filename]
markov_time_jump_string = "markov_time_jump="
markov_time_jump_label = "time_jump="

Label.string_extract(filename, markov_time_jump, start_string = markov_time_jump_string, end_string = ',')

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

# Coupling strength to Bath if exists
bath_coupling = [[] for n in filename]
bath_coupling_start = "K="
bath_coupling_label = "K="

Label.string_extract(filename, bath_coupling, start_string = bath_coupling_start,
    end_string = ',')

Label.label_build(label, '', general, end = True, concat = '')
#Label.label_build(legend, '', general, end = True, concat = '')

Label.label_build(label, length_label, length, end = True)
Label.label_build(legend, length_label, length, end = True)

Label.label_build(label, coupling_label, coupling, end = True)
Label.label_build(legend, coupling_label, coupling, end = True)

Label.label_build(label, bath_coupling_label, bath_coupling, end = True)
Label.label_build(legend, bath_coupling_label, bath_coupling, end = True)

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
#Label.label_build(legend, model_num_label, model_num, end = True)

Label.label_build(label, markov_time_jump_label, markov_time_jump, end = True)
#Label.label_build(legend, markov_time_jump_label, markov_time_jump, end = True)

Label.label_build(label, type_label, cal_type, end = True)
Label.label_build(label, task_label, cal_task, end = True)
Label.label_build(label, Total_label, Total, end = True)
Label.label_build(label, jump_label, jump, end = True)

Label.label_build(label, '', init, end = True)
#Label.label_build(legend, '', init, end = True)



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

draw1.figure_init(ymax = 1, ymin=0.001,logy = True, xmax=10)
draw1.figure_set()

#in_range = ["3.14", "1.04", "1.57", "2.09", "2.61"] # Angles
must_in_range = ["Markov Inter Random Both X Floquet","J=0.3","Model_Num=1 "]
must_not_range = []

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)



# Do linear regression
start = 3
fit_slope = []
fit_length = []
fit_slope_error = []
for i in draw1._plot_range:
    fit_x = time[i][start:]
    fit_y = [log(n) for n in spin_ave[i][start:]]
    start += 2

    slope, intercept, r_value, p_value, std_err = stats.linregress(fit_x,fit_y)

    fit_slope.append(slope)
    fit_length.append(int(length[i]))

    print '\n'
    print length_label + length[i]
    print "Fitted slope: " + str(slope)
    print "Fitted intercept: " + str(intercept)

    # Compute estimated error for slope
    mx = np.mean(fit_x)
    sx2 = np.sum((fit_x-mx)**2)
    slope_sd = std_err*sqrt(1.0/sx2)

    print "slope sd: " + str(slope_sd)
    print "R^2: " + str(r_value**2)

    fit_slope_error.append(slope_sd)

sort_index = np.argsort(fit_length)
fit_rate = [abs(n) for n in fit_slope]

draw1.errorbar(time, spin_ave, yerr = spin_err, label = legend)

pylab.legend(loc='left', ncol=1, prop={'size':14})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel(r"$\sigma^{z,l}(t)$")
pylab.xlabel("time")

#pylab.savefig("Inter_Random_Floquet_Markov_J_0_3_largest_left_spin_eigenstate_left_spin_model_num_100_markov_time_jump_10_compare_size_log.png",box_inches='tight')

pylab.figure(2)
ax=pylab.subplot(111)
ax.set_xlim(xmin=fit_length[sort_index[0]]-0.1, xmax=fit_length[sort_index[-1]]+0.1)
ax.set_ylim(ymax = max(fit_rate)+0.01, ymin = min(fit_rate)-0.01)
ax.set_yscale("log", nonposy='clip')
#pylab.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))

pylab.errorbar(fit_length, fit_rate, yerr = fit_slope_error, linewidth=0, marker='o',markersize=6)

pylab.ylabel("Decay Rate")
pylab.xlabel(r"$L$")

#pylab.savefig("Inter_Random_Floquet_Markov_J_0_3_largest_left_spin_eigenstate_left_spin_model_num_100_markov_time_jump_10_decay_rate_size_compare.pdf",box_inches='tight')

"""
min_fit_rate = 0.08
modified_fit_rate = [n - min_fit_rate for n in fit_rate]

log_fit_y = []
log_fit_x = []

for i in range(len(modified_fit_rate)):
    if fit_length[i] != 10:
        log_fit_y.append(log(modified_fit_rate[i]))
        log_fit_x.append(fit_length[i])

# Do linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(log_fit_x,log_fit_y)

print '\n'
print length_label + length[i]
print "Fitted slope for modified log fit rate: " + str(slope)
print "Fitted intercept for modified log fit rate: " + str(intercept)

# Compute estimated error for slope
mx = np.mean(fit_x)
sx2 = np.sum((fit_x-mx)**2)
slope_sd = std_err*sqrt(1.0/sx2)

print "slope sd for modified log fit rate: " + str(slope_sd)
print "R^2 for modified log fit rate: " + str(r_value**2)

from math import exp
log_fit_line = [exp(slope*n+intercept) for n in log_fit_x]



pylab.figure(3)
ax=pylab.subplot(111)
ax.set_xlim(xmin=fit_length[sort_index[0]]-0.1, xmax=fit_length[sort_index[-1]]+0.1)
#ax.set_ylim(ymax = max(fit_rate)+0.01, ymin = min(fit_rate)-0.01)
ax.set_yscale("log", nonposy='clip')
#pylab.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))

pylab.errorbar(fit_length, modified_fit_rate, yerr = fit_slope_error, linewidth=0, marker='o',markersize=6)

pylab.plot(log_fit_x, log_fit_line, color='r', linewidth=3)

pylab.ylabel("Modified Decay Rate")
pylab.xlabel(r"$L$")

#pylab.savefig("Inter_Random_Floquet_Markov_J_0_9_largest_left_spin_eigenstate_left_spin_model_num_400_log_decay_rate_subtract_0_08_size_compare.pdf",box_inches='tight')

"""

pylab.show()

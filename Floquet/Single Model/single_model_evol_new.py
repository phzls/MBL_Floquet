path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing as Draw
import label_build as Label
import numpy as np
from scipy import stats
from math import sqrt, log

# single_model file name index
name_index = 3

# long time average number
ave_num = 1000

def read_sign(filename):
    temp = np.loadtxt(filename+".txt")
    initial_spin = []
    sign = []

    for n in range(len(temp)):
        initial_spin.append(temp[n][0])
        sign.append(temp[n][1])

    return initial_spin, sign

def read_spin(filename):
    """  Read average spin values from files concerning time evolution of leftmost
         spin per mode. It is assumed there is only 1 model or 1 realization, and
         apart from the first column, which is time, each column gives averages
         under different conditions. The column with -1 time gives the index."""

    temp = np.loadtxt(filename+".txt")

    t = []
    spin_ave = []

    num_col = len(temp[0]) - 1
    if num_col > 1:
        spin_ave = [[] for i in range(num_col)]

    for n in range(len(temp)):
        if temp[n][0] < 0:
            spin_index = temp[n][1:]
        else:
            t.append(temp[n][0])
            for l in range(num_col):
                spin_ave[l].append(temp[n][l+1])
    return t, spin_ave, spin_index

filename = []
label = [] # General label for excuding when plotting

# Obtain filenames
Label.name_read("single_model_name"+str(name_index), filename, label, include_s = "Markov",
                exclude_s = "easy_full")

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
realization_string = "Run="
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

# Initial Leftmost spin z index
spin_z_index = [[] for n in filename]
spin_z_index_start = "Leftmost_spin_z_index="
spin_z_index_label = "Index="

Label.string_extract(filename, spin_z_index, start_string = spin_z_index_start, end_string = ',')

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
#Label.label_build(legend, bath_coupling_label, bath_coupling, end = True)

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

Label.label_build(label, spin_z_index_label, spin_z_index, end = True)
Label.label_build(legend, spin_z_index_label, spin_z_index, end = True)



#print label
print legend
print label

# Assume only 1 file to read
for n in range(len(filename)):
    t, spin_ave, spin_index = read_spin(filename[0])

basic_label = label[0]
basic_legend = legend[0]

time = []
for i in range(len(spin_index)):
    time.append(t)
    if i == 0:
        label[0] += " Index=" + str(int(spin_index[i]))
        legend[0] += " Index=" + str(int(spin_index[i]))
    else:
        label.append(basic_label + " Index=" + str(int(spin_index[i])))
        legend.append(basic_legend + " Index=" + str(int(spin_index[i])))

spin_ave_sub_abs = []
spin_ave_abs = []

for i in range(len(spin_ave)):
    long_time_ave = np.average(spin_ave[i][-ave_num:])
    spin_ave_sub_abs.append([abs(n - long_time_ave) for n in spin_ave[i]])
    spin_ave_abs.append([abs(n) for n in spin_ave[i]])

import pylab
draw1 = Draw.Draw()

draw1.figure_init()
draw1.figure_set()

#in_range = ["3.14", "1.04", "1.57", "2.09", "2.61"] # Angles
must_in_range = ["Markov Inter Random Both X Floquet","J=0.3"]
must_not_range = []

draw1.plot_range(label, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)


draw1.plot(time, spin_ave, label = legend)

pylab.legend(loc='upper right', ncol=2, prop={'size':12})#, bbox_to_anchor=(0.7, 1.01))
pylab.ylabel(r"$\sigma^{z,l}(t)$")
pylab.xlabel("time")

#pylab.subplots_adjust(left=0.16)

#pylab.savefig("Inter_Random_Floquet_Markov_Single_Model_J_0_3_left_spin_eigenstate_markov"
#              "_time_jump_10_compare_index_long_"+ str(ave_num)+",v" +
#              str(name_index) + ".png",box_inches='tight')

"""
# Do linear regression
start = 200000/8 + 1
end = 600000/8-1

print time[0][start], time[0][end]

fit_slope = []
fit_length = []
fit_slope_error = []
for i in draw1._plot_range:
    fit_x = time[i][start:end]
    fit_y = [log(n) for n in spin_ave_sub_abs[i][start:end]]

    slope, intercept, r_value, p_value, std_err = stats.linregress(fit_x,fit_y)

    fit_slope.append(slope)
    fit_length.append(int(length[i]))

    print '\n'
    print spin_z_index_label + spin_z_index[i]
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

# Process sign changes in evolution

f = open("single_model_name"+str(name_index)+".txt",'r')
include_s = "easy_full"
exclude_s = None
for line in f:
    include = False
    if include_s != None:
        if line.find(include_s) != -1:
            include = True
    else:
        include = True

    if exclude_s != None:
        if line.find(exclude_s) != -1:
            include = False

    if line.split()[0] != "#":
        if include:
            start = line.find(' ')
            filename_sign = line[0:start]
f.close()

ini_spin, sign = read_sign(filename_sign)

pylab.figure(3)
ax=pylab.subplot(111)
ax.set_ylim(ymin = -2, ymax = 2)
pylab.plot(ini_spin, sign, linewidth = 0, marker = '^', markersize = 10,
           markeredgecolor = 'b')
pylab.xlabel("Initial spin")

pos_change = 0
neg_change = 0
pos_num = 0
neg_num = 0
zero_pos = 0

for i in range(len(ini_spin)):
    if ini_spin[i] > 0 and sign[i] == -1:
        pos_change += 1
    elif ini_spin[i] < 0 and sign[i] == 1:
        neg_change += 1
    if sign[i] != 0:
        if ini_spin[i] > 0:
            pos_num += 1
        else:
            neg_num += 1
    else:
        zero_pos = i

print "# of pos ini which change sign: ", pos_change
print "# of neg ini which change sign: ", neg_change
print "# of pos ini: ", pos_num
print "# of neg ini: ", neg_num
print "Position of zero: ", zero_pos

pylab.savefig("Inter_Random_Floquet_Markov_Single_Model_J_0_3_left_spin_eigenstate_markov"
              "_time_jump_10_spin_evolution_sign,v" +
              str(name_index) + ".pdf",box_inches='tight')

pylab.show()



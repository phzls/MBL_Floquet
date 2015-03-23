path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing as Draw
import label_build as Label
import numpy as np

def read_spin(filename, t, spin_ave, spin_err, realization = 1):
    """  Read averages and standard deviation of entropies from files
         concerning time evolution of entropy per mode."""

    temp = np.loadtxt(filename+".txt")

    for n in range(len(temp)):
        t.append(temp[n][0])
        spin_ave.append(temp[n][1])
        spin_err.append(temp[n][2] / (realization**0.5) )

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

# Type of spin
""" I should put a more easily recordnized label in front of it """
spin_type = [[] for n in filename]
type_string = "left_spin_"
type_label = "spin_"

Label.string_extract(filename, spin_type, start_string = type_string, end_string = "_per_mode")

# Realizations
realization = [[] for n in filename]
realization_string = "Run="
realization_label = "Run="

Label.string_extract(filename, realization, start_string = realization_string, end_string = ',')

# Initial conditions
init = [[] for n in filename]
init_string = "Init_"

Label.string_extract(filename, init, start_string = init_string, end_string = ',')
Label.label_change(init)

# Total time step
Total = [[] for n in filename]
Total_string = "Total="
Total_label = "Total="

Label.string_extract(filename, Total, start_string = Total_string, end_string = ',')

# jump
jump = [[] for n in filename]
jump_string = "jump="
jump_label = "jump="

Label.string_extract(filename, jump, start_string = jump_string, end_string = ',')

# Disorder strength if exists
coupling = [[] for n in filename]
coupling_start = "lambda="
coupling_label = "J="

Label.string_extract(filename, coupling, start_string = coupling_start, end_string = ',')

# Version number if exists
version = [[] for n in filename]
version_start = ",v"
version_label = "v"

Label.string_extract(filename, version, start_string = version_start)

Label.label_build(label, '', general, end = True, concat = '')
Label.label_build(legend, '', general, end = True, concat = '')

Label.label_build(label, length_label, length, end = True)
Label.label_build(legend, length_label, length, end = True)

Label.label_build(label, coupling_label, coupling, end = True)
Label.label_build(legend, coupling_label, coupling, end = True)

Label.label_build(label, realization_label, realization, end = True)
Label.label_build(label, type_label, spin_type, end = True)
Label.label_build(label, Total_label, Total, end = True)
Label.label_build(label, jump_label, jump, end = True)
Label.label_build(label, '', init, end = True)
Label.label_build(label, version_label, version, end = True)

#print label
#print legend

class Spin_Data(object):
    def __init__(self):
        self.time = []
        self.spin_ave = []
        self.spin_err = []

spin_x = Spin_Data()
spin_y = Spin_Data()
spin_z = Spin_Data()
label_x = []
label_y = []
label_z = []
legend_x = []
legend_y = []
legend_z = []

# Dictionary for calling
spin_dict = {}
spin_dict["spin_x"] = spin_x
spin_dict["spin_y"] = spin_y
spin_dict["spin_z"] = spin_z

label_dict = {}
label_dict["spin_x"] = label_x
label_dict["spin_y"] = label_y
label_dict["spin_z"] = label_z

legend_dict = {}
legend_dict["spin_x"] = legend_x
legend_dict["spin_y"] = legend_y
legend_dict["spin_z"] = legend_z

for n in range(len(filename)):
    spin = "spin_" + spin_type[n]
    data = spin_dict[spin]
    label_dict[spin].append(label[n])
    legend_dict[spin].append(legend[n])

    data.time.append([])
    data.spin_ave.append([])
    data.spin_err.append([])

    read_spin(filename[n], data.time[-1], data.spin_ave[-1], data.spin_err[-1],float(realization[n]))

import pylab
draw1 = Draw.Draw()

draw1.figure_init(xmin = 20000, xmax = 20600)
draw1.figure_set()

s = "_x"
spin = "spin" + s
version = "v2"
run = 1
J = 0.9
L = 8

run_label = "Run=" + str(run)
if run == 1:
    run_label += " "

J_label = "J=" + str(J)
L_label = "L=" + str(L)

data = spin_dict[spin]
label_local = label_dict[spin]
legend_local = legend_dict[spin]

must_in_range = [version, run_label, J_label]
must_not_range = []

draw1.plot_range(label_local, must_in_range = must_in_range, must_not_range = must_not_range,
                 printout = True)

draw1.errorbar(data.time, data.spin_ave, yerr = data.spin_err, label = legend_local)

pylab.legend(loc='upper right', ncol=1, prop={'size':14})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel(r"$\sigma" + s + "(t)$")
pylab.xlabel("time")

J_out = list(str(J))
J_out[str(J).find('.')] = "_"
J_out = "".join(J_out)

save_name = "XXZ_Random_Floquet_" + str(L) + "_" + spin + "_J_" + J_out + "_run_" + str(run) + "_short," + version

print save_name

pylab.savefig(save_name + ".png",box_inches='tight')

pylab.show()

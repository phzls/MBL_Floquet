__author__ = 'liangshengzhang'

"""
Plot distribution of some quantities for all eigenstates
"""

path = "/Users/liangshengzhang/Google Drive/Programs/python_module/"

import sys
sys.path.append(path)

import drawing_new_1 as Draw
import label_build as Label
import numpy as np

include_s = "zz" # A string that must be included in filename construct
exclude_s = "Shift" # A string that must be excluded in filename construct
ver_num = 1

def read_file(filename):
    # Assuming different rows are for different realizations
    temp = np.loadtxt(filename+".txt")
    return temp

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

# J values
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
Label.label_build(legend, J_string, J, end = True)

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

# Find the correct version number
ver_instance = 0
for i in range(len(version)):
    if version[i] == ver_num:
        ver_instance = i
        break

print "Read file", filename[ver_instance]
dist = read_file(filename[ver_instance])

dist_max = [] # Maximum of all realizations
for n in range(len(dist)):
    dist_max.append(max(dist[n]))

instance = 0 # Which realization is plotted

import pylab

draw1 = Draw.Draw()
draw1.plot_set(logx=True,logy=True)

plot_legend = "J="+J[ver_instance] + " One Realization"

print "Average for one realization:", np.mean(dist[instance])

bin_num = 100
hist, edges = np.histogram(dist[instance], bins=bin_num)

hist_x = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]

draw1.plot(hist_x, hist, legend = plot_legend, trun=False)

pylab.legend(loc='upper right', ncol=1, prop={'size':16})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel("Frequency")
pylab.xlabel("ZZ Correlation Square")

J_val = str("J="+J[ver_instance]).replace('.','_')

save_name = "Simple_Random_Floquet_zz_corr_square_all_evec_J_" + J_val + \
            "_instance_" + str(instance)+ "_bin_num" + str(bin_num) +  "_v_" + str(ver_num)
print save_name
pylab.subplots_adjust(bottom=0.14)
pylab.savefig(save_name + "_v" + str(ver_num) + ".pdf",box_inches='tight')

draw2 = Draw.Draw()
draw2.plot_set(logx=True,logy=True)

plot_legend = "J="+J[ver_instance] + " All Realizations"

dist_all = dist.flatten()

hist, edges = np.histogram(dist_all, bins=1000)

print "Average for all realizations:", np.mean(dist_all)

hist_x = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]

draw2.plot(hist_x, hist, legend = plot_legend, trun=False)

print max(hist), draw1._ymax

pylab.legend(loc='upper right', ncol=1, prop={'size':16})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel("Frequency")
pylab.xlabel("ZZ Correlation Square")

J_val = str(J[ver_instance]).replace('.','_')
num_run = realization[ver_instance]

save_name = "Simple_Random_Floquet_zz_corr_square_all_evec_J_" + J_val + "_run_" + str(num_run)\
            + "_bin_num" + str(bin_num) + "_v_" + str(ver_num)
print save_name
pylab.subplots_adjust(bottom=0.14)
pylab.savefig(save_name + "_v" + str(ver_num) + ".pdf",box_inches='tight')


draw3 = Draw.Draw()
draw3.plot_set()

plot_legend = J[ver_instance] + " Maximunm from all Realizations"

x = range(len(dist_max))

draw3.plot(x, dist_max, legend = plot_legend, trun=False)

print max(hist), draw1._ymax

pylab.legend(loc='upper right', ncol=1, prop={'size':16})#, bbox_to_anchor=(1.1, 0.5))
pylab.ylabel("Maximum ZZ Correlation Square")
pylab.xlabel("Realization")

print "Var of max:", np.var(dist_max)

J_val = str(J[ver_instance]).replace('.','_')
num_run = realization[ver_instance]

save_name = "Simple_Random_Floquet_zz_corr_square_all_max_J_" + J_val + "_run_" + str(num_run)\
            +  "_v_" + str(ver_num)
print save_name
#pylab.subplots_adjust(bottom=0.14)
pylab.savefig(save_name + "_v" + str(ver_num) + ".pdf",box_inches='tight')

pylab.show()





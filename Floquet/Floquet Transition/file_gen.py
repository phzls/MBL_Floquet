__author__ = 'liangshengzhang'

import os
import fnmatch

# files with different J
trans_file = []

# files with same J
single_file = []

# This is not included in the file name
suffix = ".txt"

# Read all files starting with Markov
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, "*J_N=*.txt"):
#        print file[:-len(suffix)]
        trans_file.append(file[:-len(suffix)])
    if fnmatch.fnmatch(file, "*J=*.txt"):
        single_file.append(file[:-len(suffix)])

f = open("name.txt",'w')
for n in trans_file:
    print >> f, n
f.close()

f = open("single_name.txt",'w')
for n in single_file:
    print >> f, n
f.close()

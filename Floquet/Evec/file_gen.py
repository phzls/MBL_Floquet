__author__ = 'liangshengzhang'

import os
import fnmatch

trans_file = []

# This is not included in the file name
suffix = ".txt"

# Read all files starting with Markov
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, "*Floquet*.txt"):
#        print file[:-len(suffix)]
        trans_file.append(file[:-len(suffix)])

f = open("name.txt",'w')
for n in trans_file:
    print >> f, n
f.close()
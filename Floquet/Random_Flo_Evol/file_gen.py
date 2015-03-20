import os
import fnmatch

# Markov file names
markov_file = []

# This is not included in the file name
suffix = ".txt"

# Read all files starting with Markov
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, "*_spin*.txt"):
#        print file[:-len(suffix)]
        markov_file.append(file[:-len(suffix)])

f = open("spin_name.txt",'w')
for n in markov_file:
    print >> f, n
f.close()

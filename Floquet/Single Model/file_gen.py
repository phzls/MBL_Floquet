import os
import fnmatch

# single model file names
single_file = []

# Whether increase single_model_name index
name_next = True

# This is not included in the file name
txt_suffix = ".txt"
eval_phase_suffix = ",eval_phases_under_one_model,v2.txt"

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

single_num = []
# Get the current newest single_model_name file
for n in os.listdir('.'):
	if fnmatch.fnmatch(n, "single_model_name*.txt"):
		s =  find_between(n, "single_model_name", ".txt")
		single_num.append(int(s))

inc = 0
if name_next:
	inc += 1

try: 
	num = max(single_num) + inc
except ValueError:
	num = 1

print num

# Read all files starting with Markov
# Find one file with eval which generates the series of files about matrices
for n in os.listdir('.'):
    if fnmatch.fnmatch(n, "Markov*Total=120000*.txt"):
#        print n[:-len(suffix)]
        single_file.append(n[:-len(txt_suffix)])
    elif fnmatch.fnmatch(n, "*"+eval_phase_suffix):
    	single_file.append(n[:-len(eval_phase_suffix)])

f = open("single_model_name"+str(num)+".txt",'w')
for n in single_file:
    print >> f, n
f.close()

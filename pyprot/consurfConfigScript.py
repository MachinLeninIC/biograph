import sys
import json

if len(sys.argv) == 1:
    filename = 'pdbaa_list.nr'
elif len(sys.argv) == 2:
    filename = sys.argv[0]
else:
    print('too many arguments')


data={}

with open(filename) as ifile:
    for line in ifile.readlines():
        h, l = line.rstrip('\n').rstrip('.').split(':')
        l = l.replace(' ','').split(',')
        l = [x for x in l if len(x) > 0]
        l.append(h)
        data[h] = l

with open('consurf.json', 'w') as ofile:
    json.dump(data, ofile)

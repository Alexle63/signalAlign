import os
import matplotlib.pyplot as plt
import numpy
from math import *

NAME_IDX = 9
SENSITIVY_IDX = 8

path = "/home/alexma/vcf_exploration/out_results/.results"
data = []
npprecision = []
npsensitivity = []
pbprecision = []
pbsensitivity = []
indelnpegp = []
indelnpegs = []
indelpbegp = []
indelpbegs = []
indelnpsgp = []
indelnpsgs = []
indelpbsgp = []
indelpbsgs = []
indelnpngp = []
indelnpngs = []
indelpbngp = []
indelpbngs = []
npegs = []
npegp = []
npsgs = []
npsgp = []
npngs = []
npngp = []
pbegs = []
pbegp = []
pbsgs = []
pbsgp = []
pbngp = []
pbngs = []
size = []
names = []
nameholder = 0
ballsizecounter = 0
counter2 = 0
counter3 = 0
counter4 = 0
position = 0
largestsize = 0
pbsize = []
npsize = []
file_name = os.path.join(path, "summary.no_gap_deletion.txt")
my_file = open(file_name)
for line in my_file:
    line=line.replace("NaN", "0.0")
    data.append(line.strip().split('\t'))
    print data
for d in data:
    name = d[NAME_IDX]
    names.append(name.split('-'))
    print(names)
for x in data:
    if nameholder == nameholder:
        size.append(int(x[2]))
        if names[nameholder][0] == "np":
            npsensitivity.append(float(x[SENSITIVY_IDX]))
            npprecision.append(float(x[6]))
            if names[nameholder][1] == "eg":
                npegs.append(float(x[SENSITIVY_IDX]))
                npegp.append(float(x[6]))
                if nameholder % 2 is 0:
                    indelnpegp.append(float(x[6]))
                    indelnpegs.append(float(x[SENSITIVY_IDX]))
            if names[nameholder][1] == "sg":
                npsgs.append(float(x[SENSITIVY_IDX]))
                npsgp.append(float(x[6]))
                if nameholder % 2 is 0:
                    indelnpsgp.append(float(x[6]))
                    indelnpsgs.append(float(x[SENSITIVY_IDX]))
            if names[nameholder][1] == 'ng':
                npngs.append(float(x[SENSITIVY_IDX]))
                npngp.append(float(x[6]))
                if nameholder % 2 is 0:
                    indelnpngp.append(float(x[6]))
                    indelnpngs.append(float(x[SENSITIVY_IDX]))

        else:
            pbsensitivity.append(float(x[SENSITIVY_IDX]))
            pbprecision.append(float(x[6]))
            if names[nameholder][1] == "eg":
                pbegs.append(float(x[SENSITIVY_IDX]))
                pbegp.append(float(x[6]))
                if nameholder % 2 is 0:
                    indelpbegp.append(float(x[6]))
                    indelpbegs.append(float(x[SENSITIVY_IDX]))
            if names[nameholder][1] == "sg":
                pbsgs.append(float(x[SENSITIVY_IDX]))
                pbsgp.append(float(x[6]))
                if nameholder % 2 is 0:
                    indelpbsgp.append(float(x[6]))
                    indelpbsgs.append(float(x[SENSITIVY_IDX]))
            if names[nameholder][1] == "ng":
                pbngs.append(float(x[SENSITIVY_IDX]))
                pbngp.append(float(x[6]))
                if nameholder % 2 is 0:
                    indelpbngp.append(float(x[6]))
                    indelpbngs.append(float(x[SENSITIVY_IDX]))
    nameholder+=1
for y in size:
    if size[ballsizecounter] > largestsize:
        largestsize = size[ballsizecounter]
    ballsizecounter+=1
for z in size:
    size[counter2] = size[counter2]*500
    size[counter2] = size[counter2]/largestsize
    counter2+=1
for s in names:
    if counter3 == counter3:
        if names[counter3][0] == "np":
            npsize.append(size[counter4])
        else:
            pbsize.append(size[counter4])
        counter4+=1
    counter3+=1
plt.subplot(2,1,1)
plt.xlabel("Precision")
plt.ylabel("Sensitivity")
plt.title("Nanopore")
plt.scatter(npegp, npegs, s = 400, c = "b", marker = "*", label = '\nE gap\n', alpha = 0.5)
plt.scatter(indelnpegp, indelnpegs, s= 450, edgecolors= 'r', marker='*', facecolor = 'None')
plt.scatter(npsgp,npsgs, s = 400, c = "b", marker = "o",label = '\nS gap\n', alpha = 0.5)
plt.scatter(indelnpsgp, indelnpsgs, s= 450, edgecolors= 'r', marker='o', facecolor = 'None')
plt.scatter(npngp, npngs, s = 400, c= "b", marker = "^",label = '\nN gap\n', alpha = 0.5)
plt.scatter(indelnpngp, indelnpngs, s= 450, edgecolors= 'r', marker='^', facecolor = 'None')
axes = plt.gca()
axes.set_xlim(0.45,1)
axes.set_ylim(0.45,1)
plt.plot(npegp,npegs,'b--')
plt.plot(npsgp, npsgs, 'b--')
plt.plot(npngp, npngs, 'b--')
plt.plot(0.55,0.65, 'r-', label = 'Indels')
plt.subplots_adjust(left=0.125, bottom=0.05, right=0.9, top=0.96,
                wspace=0.2, hspace=0.2)
plt.legend(loc = 2, markerscale = 1, ncol=2)
plt.subplot(212)
plt.scatter(pbegp, pbegs, s = 400, c = "g", marker = "*", label = "\nE gap\n", alpha = 0.5)
plt.scatter(indelpbegp, indelpbegs, s= 450, edgecolors= 'r', marker='*', facecolor = 'None')
plt.scatter(pbsgp,pbsgs, s = 400, c = "g", marker = "o", label = "\nS gap\n", alpha = 0.5)
plt.scatter(indelpbsgp, indelpbsgs, s= 450, edgecolors= 'r', marker='o', facecolor = 'None')
plt.scatter(pbngp, pbngs, s = 400, c= "g", marker = "^", label = "\nN gap\n", alpha = 0.5)
plt.scatter(indelpbngp, indelpbngs, s= 450, edgecolors= 'r', marker='^', facecolor = 'None')
axes = plt.gca()
axes.set_xlim(0.45,1)
axes.set_ylim(0.45,1)
plt.plot(pbegp,pbegs,'g--')
plt.plot(pbsgp, pbsgs, 'g--')
plt.plot(pbngp, pbngs, 'g--')
plt.plot(0.97,0.9, 'r-', label = 'Indels')
plt.xlabel("Precision")
plt.ylabel("Sensitivity")
plt.title("Pacbio")
plt.legend(loc = 2, markerscale = 1, ncol = 2)
plt.show()


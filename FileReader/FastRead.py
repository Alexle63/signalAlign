import os
import matplotlib.pyplot as plt
import numpy as np
from math import *

FILE_DATA = []
CHROMOSOME_IDX = 1
CHROMOSOME = []
CHROMOSOME_VALUE = []
CHROMOSOME_COUNT = []
counter = 0
place2 = 0
counter2= 0
PLACE_COUNTER = 0
place = 0
path = "/home/alexma/dev2/signalAlign/output/tempFiles_errorCorrection/step_0/"
file_name = os.path.join(path, "4f614987-de21-472a-b58e-71faff4bb369.sm.forward.tsv")
my_file = open(file_name)
for line in my_file:
    FILE_DATA.append(line.strip().split('\t'))
CHROMOSOME_VALUE.append(FILE_DATA[0][1])
for x in FILE_DATA:
    CHROMOSOME.append(x[CHROMOSOME_IDX])
    if CHROMOSOME_VALUE[place] != CHROMOSOME[counter]:
        CHROMOSOME_VALUE.append(CHROMOSOME[counter])
        place+=1
        counter+=1
    else:
         counter+=1
CHROMOSOME_COUNT.append(0)
print(CHROMOSOME)
for y in CHROMOSOME:
    if CHROMOSOME_VALUE[place2] != CHROMOSOME[counter2]:
        CHROMOSOME_COUNT.append(1)
        counter2+=1
        place2+=1
    else:
        CHROMOSOME_COUNT[place2]+=1
        counter2+=1
PLACE_COUNTER = range(len(CHROMOSOME_COUNT))
#plt.plot(CHROMOSOME_COUNT)
#plt.bar(PLACE_COUNTER,CHROMOSOME_COUNT)
#plt.xticks(PLACE_COUNTER, CHROMOSOME_VALUE)
#plt.show()

my_data = dict()
with open(file_name) as my_file:
    for line in my_file:
        line = line.split("\t")
        index = int(line[1])
        if index not in my_data:
            my_data[index] = list()
        my_data[index].append(float(line[12]))
print(my_data)
x_coords = list(my_data.keys())
x_coords.sort()
print(x_coords)
y_coords = list()
print(y_coords)
for x in x_coords:
    y_coords.append(len(my_data[x]))
print(y_coords)

PLACE_COUNTER = range(len(y_coords))
plt.bar(PLACE_COUNTER, y_coords)
plt.xticks(PLACE_COUNTER, x_coords)
plt.show()

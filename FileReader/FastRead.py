import os
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser
import argparse
from math import *

FILE_DATA = []
CHROMOSOME_IDX = 1
CHROMOSOME = []
CHROMOSOME_VALUE = []
CHROMOSOME_COUNT = []
MINIMUM_X = []
MINIMUM_Y = []
MIN_EVENT_X = []
MIN_EVENT_Y = []
avg_probs = []
PRINT_X_ARRAY = []
prob = []
PROB_LIST = []
steps = []
CHROM_POSITION = []
numindex2 = []
#plt.plot(CHROMOSOME_COUNT)
#plt.bar(PLACE_COUNTER,CHROMOSOME_COUNT)
#plt.xticks(PLACE_COUNTER, CHROMOSOME_VALUE)
#plt.show()

def parse_args(args=None):
     parser = argparse.ArgumentParser(description=__doc__)
     parser.add_argument('--start_position', '-s', action='store',
                         required=False, type=int, default=0,
                         help="start index of reference points/event positions")
     parser.add_argument('--range', '-r', action='store',
                           required=True, type=int, default=250,
                          help="number of reference points/event positions in 1 step")
     parser.add_argument('--tsv_file', '-f', action='store',
                         required=True, type=str, default=None,
                         help="TSV file path to be graphed")
     parser.add_argument('--position', '-p', action='store',
                         required=False, type=int, default=1,
                         help="index for which column for x axis. 1 = Reference positions, 5 = Event positions")
     args = parser.parse_args()
     return args


def main():
    args = parse_args()
    startpos = args.start_position
    rang = args.range
    tsvfile = args.tsv_file
    pos = args.position
    CHROMOSOME_POSITION_COUNTER = 0
    INDEX = pos
    average = 0
    PRINT_COUNTER = 0
    switcher = 0
    pertick = int((rang+50)/100)
    length = 0
    CURRENT_PROBABILITY = 0
    POSITION_COUNT = 0
    place2 = 0
    value_counter = 0
    PLACE_COUNTER = 0
    counter3 = startpos  # change values
    maxima = rang+counter3  # change values
    BOT_START = 0  # change values
    BOT_END = rang  # change values
    place = 0
    print(maxima)
    my_data = dict()
    with open(tsvfile) as my_file:
        for line in my_file:
            line = line.split("\t")
            index = int(line[INDEX])
            if index not in my_data:
                my_data[index] = list()
            my_data[index].append(float(line[12]))
    x_coords = list(my_data.keys())
    x_coords.sort()
    y_coords = list()
    for x in x_coords:
        y_coords.append(len(my_data[x]))
    while counter3 < maxima:
        MINIMUM_X.append(x_coords[counter3])
        MINIMUM_Y.append(y_coords[counter3])
        counter3+=1
    for x in my_data:
        PROB_LIST.append(my_data[x_coords[CHROMOSOME_POSITION_COUNTER]])
        CHROMOSOME_POSITION_COUNTER+=1
    #PLACE_COUNTER = range(len(y_coords))
    #plt.bar(PLACE_COUNTER, y_coords)
    #print(PLACE_COUNTER)
    #plt.xticks(PLACE_COUNTER, x_coords)
    #plt.show()
    counter11 = startpos
    counter12 = 0
    maximus = []
    for w in MINIMUM_X:
        maximus.append(max(list(my_data.values())[counter11]))
        counter11+=1
        maximus[counter12] = 100*(1-maximus[counter12])
        counter12+=1
    counter3 = startpos
    counter4 = 0
    counter5 = 0
    total = []
    top = 0
    final = 0
    first = []
    for z in MINIMUM_X:
        average = sum(PROB_LIST[counter3])
        length = len(PROB_LIST[counter3])
        average = average/length
        while value_counter <= ((MINIMUM_Y[counter4])-1):
            prob.append(PROB_LIST[counter3][value_counter])
            value_counter+=1
        while counter5 <=((MINIMUM_Y[counter4])-1):
            total.append(abs(average - prob[counter5]))
            counter5+=1
        top = sum(total)
        final = (top/MINIMUM_Y[counter4])
        final = 2*final
        first.append(final)
        total.clear()
        prob.clear()
        counter5 = 0
        counter4+=1
        value_counter = 0
        counter3+=1
    value_counter2 = 0
    probabilities = []
    counter6 = startpos
    counter7 = 0
    counter10 = 0
    reciprocals = []
    for a in MINIMUM_X:
        while value_counter2<=((MINIMUM_Y[counter7])-1):
            probabilities.append(PROB_LIST[counter6][value_counter2])
            value_counter2+=1
            reciprocals.append(1/probabilities[counter10])
            counter10+=1
        value_counter2 = 0
        counter6+=1
        counter7+=1
    value_counter3 = 0
    counter8 = 0
    counter9 = 0
    bottom = []
    final2 = []
    hentropy = []
    for s in MINIMUM_X:
        num = 0
        num2 = 0
        while value_counter3<=((MINIMUM_Y[counter8])-1):
            num = reciprocals[counter9]
            num2+=num
            value_counter3+=1
            counter9+=1
        bottom.append(num2)
        final2.append(1-(MINIMUM_Y[counter8]/bottom[counter8]))
        value_counter3 = 0
        counter8+=1
    value_counter4 = 0
    counter14 = 0
    counter13 = 0
    for y in MINIMUM_X:
        num3 = 0
        num4 = 0
        while value_counter4<=((MINIMUM_Y[counter13])-1):
            num3 = probabilities[counter14]*(np.log2(probabilities[counter14]))
            num4+=num3
            counter14+=1
            value_counter4+=1
        hentropy.append(num4*-1)
        counter13+=1
        value_counter4 = 0
    while PRINT_COUNTER< len(MINIMUM_X):
        PRINT_X_ARRAY.append(MINIMUM_X[PRINT_COUNTER])
        PRINT_COUNTER+=pertick#change value
    my_file.close()
    f, (sub1, sub2, sub3, sub4, sub5) = plt.subplots(5,1, sharex= True)
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub1.bar(PLACE_COUNTER, first)
    sub1.set(ylabel = 'Pesout Theorem\n-------------------\n')
    plt.subplots_adjust(left=0.06, bottom=0.12, right=1, top=.98,
                    wspace=0.2, hspace=0.05)
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub2.bar(PLACE_COUNTER, final2)
    sub2.set(ylabel = 'Harmonic Mean\n-------------------\n')
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub3.bar(PLACE_COUNTER, MINIMUM_Y)
    sub3.set(ylabel = 'Position Count\n-------------------\n\n')
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub4.bar(PLACE_COUNTER, maximus)
    sub4.set(ylabel = 'Percent Uncertainty\n-------------------\n')
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub5.bar(PLACE_COUNTER, hentropy)
    sub5.set(ylabel = 'Huffman Entropy\n-------------------\n')
    plt.xticks(np.arange(BOT_START,BOT_END,pertick),PRINT_X_ARRAY, rotation= 90)#change value
    plt.xlabel("Chromosome Position")
    plt.show()


if __name__ == "__main__":
    main()

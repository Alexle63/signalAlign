import matplotlib.pyplot as plt
import numpy as np
import argparse
import itertools
import copy
import glob
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
SHORT_PROBABILITIES = []
NUCLEOTIDES = []
highestprobarray = []
totalofaverage = []
pesouttheorem = []
probabilities = []
reciprocals = []
harmonicbottom = []
harmonicmeanfinal = []
hentropy = []
Alist = []
Clist = []
Glist = []
Tlist = []
Atotal = []
Ctotal = []
Gtotal = []
Ttotal = []
chromosomeholder = []
chromosomepos = dict()
chromosomeprob = dict()
def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start_position', '-s', action='store',
                        required=False, type=int, default=0,
                        help="start index of reference points/event positions")
    parser.add_argument('--range', '-r', action='store',
                        required=True, type=int, default=250,
                        help="number of reference points/event positions in 1 step")
    parser.add_argument('--position', '-p', action='store',
                        required=False, type=str, default="reference",
                        help="index for which column for x axis. 1 = Reference positions, 5 = Event positions")
    parser.add_argument('--output', '-o', action='store', required=False, type=str, default=None,
                        help="Path to output destination")
    parser.add_argument('--input', '-i', action='store', required=True, type=str, default=None,
                        help="Path to input file")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    tsvFilePaths = glob.glob(args.input + '/*.tsv')
    print(tsvFilePaths)
    startpos = args.start_position
    rang = args.range
    pos = args.position
    opath = args.output + ".png"
    fileIndex = 0
    for x in tsvFilePaths:
        tsvfile = tsvFilePaths[fileIndex]
        startup(startpos, rang, pos, opath, tsvfile)
        fileIndex += 1
    #############################################################################

def startup(startpos, rang, pos, opath, tsvfile):
    if pos == "reference":
        INDEX = 1
    if pos == "event":
        INDEX = 5
    betweenseq = ""
    nucleocounter = 0
    pertick = int(rang / 100)
    counter = 0
    counter3 = 0
    maxima = rang + counter3
    BOT_START = 0
    BOT_END = rang
    sub_string = 0
    marker_index = 0
    sub_counter = 0
    my_data = dict()
    my_data2 = dict()
    my_data3 = dict()
    my_data4 = dict()
    huffmandict = dict()
    huffmandict2 = dict()
    with open(tsvfile) as my_file:
        for line in my_file:
            line = line.split("\t")
            index = int(line[INDEX])
            sub_string = len(line[15])-1
            NUCLEOTIDES.append(line[15][:sub_string])
            if index not in my_data:
                my_data[index] = list()
                my_data2[index] = list()
                my_data3[index] = list()
                my_data4[index] = list()
                chromosomepos[index] = list()
                huffmandict[index] = list()
                huffmandict2[index] = list()
            betweenseq+=str(line[12])
            betweenseq += "-"
            betweenseq += NUCLEOTIDES[nucleocounter]
            nucleocounter+=1
            my_data4[index].append(betweenseq)
            betweenseq = ""
            my_data[index].append(float(line[12]))
    my_file.close()
    x_coords = list(my_data.keys())
    x_coords.sort()
    for a in range(len(my_data)):
        endvalue = len(my_data) - a
        if endvalue < len(NUCLEOTIDES[a]):
            limit = len(NUCLEOTIDES[a]) - endvalue
            while counter < len(NUCLEOTIDES[a]) - limit:
                while sub_counter < len(my_data4[x_coords[a]]):
                    marker_index = my_data4[x_coords[a]][sub_counter].index('-')
                    my_data2[x_coords[a + counter]].append(str(my_data4[x_coords[a]][sub_counter][:marker_index] + "-" + my_data4[x_coords[a]][sub_counter][counter+1+marker_index:]))
                    sub_counter+=1
                sub_counter = 0
                counter += 1
        else:
            while counter < len(NUCLEOTIDES[a]):
                while sub_counter < len(my_data4[x_coords[a]]):
                    marker_index = my_data4[x_coords[a]][sub_counter].index('-')
                    my_data2[x_coords[a + counter]].append(str(my_data4[x_coords[a]][sub_counter][:marker_index] + "-" + my_data4[x_coords[a]][sub_counter][counter+1+marker_index:]))
                    sub_counter+=1
                sub_counter = 0
                counter += 1
        counter = 0
    #for v in range(len(my_data2)):
    #    SHORT_PROBABILITIES.append(list(itertools.chain.from_iterable(list(my_data2.values())[v])))
    #    for a in range(len(SHORT_PROBABILITIES[v])):
    #        my_data3[x_coords[v]].append(SHORT_PROBABILITIES[v][a])
    x_coords = list(my_data2.keys())
    x_coords.sort()
    y_coords = list()
    for x in x_coords:
        y_coords.append(len(my_data2[x]))

    SHORT_PROBABILITIES.clear()
    counter = 0
    for a in range(len(my_data)):
        endvalue = len(my_data) - a
        if endvalue < len(NUCLEOTIDES[a]):
            limit = len(NUCLEOTIDES[a]) - endvalue
            while counter < len(NUCLEOTIDES[a]) - limit:
                huffmandict[x_coords[a + counter]].append(my_data[x_coords[a]])
                counter += 1
        else:
            while counter < len(NUCLEOTIDES[a]):
                huffmandict[x_coords[a + counter]].append(my_data[x_coords[a]])
                counter += 1
        counter = 0
    for v in range(len(huffmandict)):
        SHORT_PROBABILITIES.append(list(itertools.chain.from_iterable(list(huffmandict.values())[v])))
        for a in range(len(SHORT_PROBABILITIES[v])):
            huffmandict2[x_coords[v]].append(SHORT_PROBABILITIES[v][a])
    x_coords2 = list(huffmandict2.keys())
    x_coords2.sort()
    y_coords2 = list()
    for x in x_coords2:
        y_coords2.append(len(huffmandict2[x]))
    maxima = len(x_coords2)
    while counter3 < maxima:
        MINIMUM_X.append(x_coords2[counter3])
        MINIMUM_Y.append(y_coords2[counter3])
        counter3 += 1
    chromosomesort(my_data2, 0, x_coords, 0, 0, 0, 0, 0, 0, 0, 0)
    highestprobability(0, 0, huffmandict2)
    sumofaverage(huffmandict2, x_coords2, 0, 0, 0, 0, 0)
    harmonicmean(0, 0, 0, 0, 0, 0, 0)
    huffmanentropy(0, 0, 0)
    printcalculations(pertick, 0)
    #printplot(BOT_START, BOT_END, pertick, opath)
    for printtrack in range(len(x_coords)):
        print("%-10d \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \t" % (x_coords[printtrack], hentropy[printtrack], Alist[printtrack], Clist[printtrack], Glist[printtrack], Tlist[printtrack]))
    FILE_DATA.clear()
    CHROMOSOME.clear()
    CHROMOSOME_VALUE.clear()
    CHROMOSOME_COUNT.clear()
    MINIMUM_X.clear()
    MINIMUM_Y.clear()
    MIN_EVENT_X.clear()
    MIN_EVENT_Y.clear()
    avg_probs.clear()
    PRINT_X_ARRAY.clear()
    prob.clear()
    PROB_LIST.clear()
    steps.clear()
    CHROM_POSITION.clear()
    numindex2.clear()
    SHORT_PROBABILITIES.clear()
    NUCLEOTIDES.clear()
    highestprobarray.clear()
    totalofaverage.clear()
    pesouttheorem.clear()
    probabilities.clear()
    reciprocals.clear()
    harmonicbottom.clear()
    harmonicmeanfinal.clear()
    hentropy.clear()
    Alist.clear()
    Clist.clear()
    Glist.clear()
    Tlist.clear()
    Atotal.clear()
    Ctotal.clear()
    Gtotal.clear()
    Ttotal.clear()
    chromosomeholder.clear()
    chromosomepos.clear()
    chromosomeprob.clear()

    #############################################################################

def chromosomesort(my_data2, substringcounter, x_coords, subcounter, basecounter, substringposcounter, Atotal, Ctotal, Gtotal, Ttotal, wholetotal):
    chromosomeprob['A'] = list()
    chromosomeprob['C'] = list()
    chromosomeprob['G'] = list()
    chromosomeprob['T'] = list()
    for r in x_coords:
        while subcounter < len(my_data2[x_coords[substringcounter]]):
            substringpos = my_data2[x_coords[substringcounter]][subcounter].index('-')
            currentbase = my_data2[x_coords[substringcounter]][subcounter][1+substringpos:basecounter+2+substringpos]
            if currentbase == "A":
                chromosomeprob['A'].append(float(my_data2[x_coords[substringcounter]][subcounter][:substringpos]))
                chromosomeprob['C'].append(0)
                chromosomeprob['G'].append(0)
                chromosomeprob['T'].append(0)
            if currentbase == "C":
                chromosomeprob['A'].append(0)
                chromosomeprob['C'].append(float(my_data2[x_coords[substringcounter]][subcounter][:substringpos]))
                chromosomeprob['G'].append(0)
                chromosomeprob['T'].append(0)
            if currentbase == "G":
                chromosomeprob['A'].append(0)
                chromosomeprob['C'].append(0)
                chromosomeprob['G'].append(float(my_data2[x_coords[substringcounter]][subcounter][:substringpos]))
                chromosomeprob['T'].append(0)
            if currentbase == "T":
                chromosomeprob['A'].append(0)
                chromosomeprob['C'].append(0)
                chromosomeprob['G'].append(0)
                chromosomeprob['T'].append(float(my_data2[x_coords[substringcounter]][subcounter][:substringpos]))
            subcounter+=1
        chromosomepos[x_coords[substringcounter]] = copy.deepcopy(chromosomeprob)
        chromosomeprob['A'].clear()
        chromosomeprob['G'].clear()
        chromosomeprob['C'].clear()
        chromosomeprob['T'].clear()
        subcounter = 0
        substringcounter+=1
    while substringposcounter < len(x_coords):
        Atotal = sum(chromosomepos[x_coords[substringposcounter]]['A'])
        Ctotal = sum(chromosomepos[x_coords[substringposcounter]]['C'])
        Gtotal = sum(chromosomepos[x_coords[substringposcounter]]['G'])
        Ttotal = sum(chromosomepos[x_coords[substringposcounter]]['T'])
        wholetotal = Atotal + Ctotal + Gtotal + Ttotal
        Alist.append(Atotal/wholetotal)
        Clist.append(Ctotal/wholetotal)
        Glist.append(Gtotal/wholetotal)
        Tlist.append(Ttotal/wholetotal)
        substringposcounter+=1

def highestprobability(startpositioncounter, highestprobabilitycounter, probabilitydict):
    for w in MINIMUM_X:
        highestprobarray.append(max(list(probabilitydict.values())[startpositioncounter]))
        startpositioncounter += 1
        highestprobarray[highestprobabilitycounter] = 100 * (1 - highestprobarray[highestprobabilitycounter])
        highestprobabilitycounter += 1
    #########################################################################

def sumofaverage(probabilitydict2, x_coords, chromosomepositioncounter, startpositioncounter2, averagearraycounter,
                 numinaveragecounter, value_counter):
    for x in range(len(probabilitydict2)):
        PROB_LIST.append(probabilitydict2[x_coords[chromosomepositioncounter]])
        chromosomepositioncounter += 1
    for z in range(len(MINIMUM_X)):
        average = sum(PROB_LIST[startpositioncounter2])
        length = len(PROB_LIST[startpositioncounter2])
        average = average / length
        while value_counter <= ((MINIMUM_Y[averagearraycounter]) - 1):
            prob.append(PROB_LIST[startpositioncounter2][value_counter])
            value_counter += 1
        while numinaveragecounter <= ((MINIMUM_Y[averagearraycounter]) - 1):
            totalofaverage.append(abs(average - prob[numinaveragecounter]))
            numinaveragecounter += 1
        top = sum(totalofaverage)
        final = (top / MINIMUM_Y[averagearraycounter])
        final = 2 * final
        pesouttheorem.append(final)
        totalofaverage.clear()
        prob.clear()
        numinaveragecounter = 0
        averagearraycounter += 1
        value_counter = 0
        startpositioncounter2 += 1
    #############################################################################

def harmonicmean(startpositioncounter3, value_counter2, value_counter3, positionprobcounter, positionprobcounter2,
                 harmonicmeancounter, harmonicreciprocalcounter):
    for a in MINIMUM_X:
        while value_counter2 <= ((MINIMUM_Y[positionprobcounter]) - 1):
            probabilities.append(PROB_LIST[startpositioncounter3][value_counter2])
            value_counter2 += 1
            reciprocals.append(1 / probabilities[harmonicmeancounter])
            harmonicmeancounter += 1
        value_counter2 = 0
        startpositioncounter3 += 1
        positionprobcounter += 1
    for s in MINIMUM_X:
        totalreciprocal = 0
        while value_counter3 <= ((MINIMUM_Y[positionprobcounter2]) - 1):
            currentreciprocal = reciprocals[harmonicreciprocalcounter]
            totalreciprocal += currentreciprocal
            value_counter3 += 1
            harmonicreciprocalcounter += 1
        harmonicbottom.append(totalreciprocal)
        harmonicmeanfinal.append(1 - (MINIMUM_Y[positionprobcounter2] / harmonicbottom[positionprobcounter2]))
        value_counter3 = 0
        positionprobcounter2 += 1

    ####################################################################################

def huffmanentropy(value_counter4, huffmanprobcounter, huffmanentropycounter):
    for y in MINIMUM_X:
        num4 = 0
        while value_counter4 <= ((MINIMUM_Y[huffmanprobcounter]) - 1):
            num3 = probabilities[huffmanentropycounter] * (np.log2(probabilities[huffmanentropycounter]))
            num4 += num3
            huffmanentropycounter += 1
            value_counter4 += 1
        if num4 == 0:
            hentropy.append(0)
        else:
            hentropy.append((num4 * -1)/10)
        huffmanprobcounter += 1
        value_counter4 = 0
    ###############################################################################

def printcalculations(pertick, PRINT_COUNTER):
    while PRINT_COUNTER < len(MINIMUM_X):
        PRINT_X_ARRAY.append(MINIMUM_X[PRINT_COUNTER])
        PRINT_COUNTER += pertick

def printplot(BOT_START, BOT_END, pertick, opath):
    ############################################################################
    f, (sub1, sub2, sub3, sub4, sub5) = plt.subplots(5, 1, sharex=True)
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub1.bar(PLACE_COUNTER, pesouttheorem)
    sub1.set(ylabel='Pesout Theorem\n-------------------\n')
    plt.subplots_adjust(left=0.06, bottom=0.12, right=1, top=.98,
                        wspace=0.2, hspace=0.05)
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub2.bar(PLACE_COUNTER, harmonicmeanfinal)
    sub2.set(ylabel='Harmonic Mean\n-------------------\n')
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub3.bar(PLACE_COUNTER, MINIMUM_Y)
    sub3.set(ylabel='Position Count\n-------------------\n\n')
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub4.bar(PLACE_COUNTER, highestprobarray)
    sub4.set(ylabel='Percent Uncertainty\n-------------------\n')
    PLACE_COUNTER = range(len(MINIMUM_Y))
    sub5.bar(PLACE_COUNTER, hentropy)
    sub5.set(ylabel='Huffman Entropy\n-------------------\n')
    plt.xticks(np.arange(BOT_START, BOT_END, pertick), PRINT_X_ARRAY, rotation=90)  # change value
    plt.xlabel("Position")
    plt.show()
    plt.savefig(opath)

if __name__ == "__main__":
    main()


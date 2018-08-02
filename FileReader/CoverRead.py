
import matplotlib.pyplot as plt
import numpy as np
import argparse
import itertools

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
                       required=False, type=str, default="reference",
                       help="index for which column for x axis. 1 = Reference positions, 5 = Event positions")
   parser.add_argument('--output', '-o', action='store', required=True, type=str, default=None,
                       help="Path to output destination")
   args = parser.parse_args()
   return args


def main():
   args = parse_args()
   startpos = args.start_position
   rang = args.range
   tsvfile = args.tsv_file
   pos = args.position
   opath = args.output + ".png"
   if pos == "reference":
       INDEX = 1
   if pos == "event":
       INDEX = 5
   pertick = int(rang / 100)
   counter = 0
   counter3 = startpos
   maxima = rang + counter3
   BOT_START = 0
   BOT_END = rang
   my_data = dict()
   my_data2 = dict()
   my_data3 = dict()
   with open(tsvfile) as my_file:
       for line in my_file:
           line = line.split("\t")
           index = int(line[INDEX])
           NUCLEOTIDES.append(line[9])
           if index not in my_data:
               my_data[index] = list()
               my_data2[index] = list()
               my_data3[index] = list()
           my_data[index].append(float(line[12]))
   my_file.close()
   x_coords = list(my_data.keys())
   x_coords.sort()
   for a in range(len(my_data)):
       endvalue = len(my_data) - a
       if endvalue < len(NUCLEOTIDES[a]):
           limit = len(NUCLEOTIDES[a]) - endvalue
           while counter < len(NUCLEOTIDES[a]) - limit:
               my_data2[x_coords[a + counter]].append(my_data[x_coords[a]])
               counter += 1
       else:
           while counter < len(NUCLEOTIDES[a]):
               my_data2[x_coords[a + counter]].append(my_data[x_coords[a]])
               counter += 1
       counter = 0
   for v in range(len(my_data2)):
       SHORT_PROBABILITIES.append(list(itertools.chain.from_iterable(list(my_data2.values())[v])))
       for a in range(len(SHORT_PROBABILITIES[v])):
           my_data3[x_coords[v]].append(SHORT_PROBABILITIES[v][a])
   x_coords = list(my_data3.keys())
   x_coords.sort()
   y_coords = list()
   print(my_data3)
   for x in x_coords:
       y_coords.append(len(my_data3[x]))
   while counter3 < maxima:
       MINIMUM_X.append(x_coords[counter3])
       MINIMUM_Y.append(y_coords[counter3])
       counter3 += 1

   highestprobability(startpos, 0, my_data3)
   sumofaverage(my_data3, x_coords, 0, startpos, 0, 0, 0)
   harmonicmean(startpos, 0, 0, 0, 0, 0, 0)
   huffmanentropy(0, 0, 0)
   printcalculations(pertick, 0)
   printplot(BOT_START, BOT_END, pertick, opath)
   #############################################################################


def highestprobability(startpositioncounter, highestprobabilitycounter, probabilitydict):
   for w in MINIMUM_X:
       highestprobarray.append(max(list(probabilitydict.values())[startpositioncounter]))
       startpositioncounter += 1
       highestprobarray[highestprobabilitycounter] = 100 * (1 - highestprobarray[highestprobabilitycounter])
       highestprobabilitycounter += 1
   #########################################################################


def sumofaverage(probabilitydict2, x_coords, chromosomepositioncounter, startpositioncounter2, averagearraycounter,
                numinaveragecounter, value_counter):
   for x in probabilitydict2:
       PROB_LIST.append(probabilitydict2[x_coords[chromosomepositioncounter]])
       chromosomepositioncounter += 1
   for z in MINIMUM_X:
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

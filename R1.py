#!/usr/bin/env python

#Converting characters into a score
def convert_phred(letter):
    # YOUR CODE HERE
    phred_score = ord(letter) -33  #ord is a function that turns letters into ASCII 
    return(phred_score)

import gzip 

f_1 = '/projects/bgmp/dchilin/demultiplexing/1294_S1_L008_R1_001.fastq.gz'
def populate_array(f_1):
    LN = 0
    mean_scores = [0.0] * 101
    with gzip.open(f_1, 'rt') as fh:
        i = 0
        for line in fh:
            i+=1
            line = line.strip('\n')
            if i % 4 == 0: #this gives us every 4th line which contains the ASCII letters for the phred score
            #everything above this is to opening a file, reading one line at at time, and grabbing the 4th line   
                index = 0 #index in the array, it should start at 0.
                for character in line: #for every ASCII character in the the 4th line
                    mean_scores[index] += (convert_phred(character)) #convert each character to a phred score number
                    #converted scores will be incremented into the mean score array.
                    index += 1   #idex will increment onto the next index within the array so that new character can be converted 
                #for b, x in enumerate(line): #for every letter(x) in the line
                    #mean_scores[b] += (convert_phred(x)) #convert each ASCII letter to a phred score number           
            LN += 1
    return mean_scores, LN

#make sure file is in the directory you are working on. That includes jupyter notebook
mean_scores, LN = populate_array(f_1)         #a is mean_scores b is LN       
#print(mean_scores)    

mean_scores, LN = populate_array(f_1)
print("QS Distribution per nulceotide-Read 1")
print("# Base Pair\tMean Quality Score")
for index, scores in enumerate(mean_scores): #for index in the mean score array
    #enumerate will number out  the score from the array
    total = (LN/4)
    ms_div_total = mean_scores[index]/total
    mean_scores[index] = float(ms_div_total)
    print(str(index)+"\t"+str(mean_scores[index]))
    




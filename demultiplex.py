#!/usr/bin/env python

#SBATCH --partition=short    ### Partition
#SBATCH --job-name=demult  ### Job Name
#SBATCH --output=demult.out    ### File in which to store job output
#SBATCH --time=1-00:00:00   ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1           ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=14 ### Number of tasks to be launcged per Node


#################################################################

#create a new working directory for this assignment on your local computer so you can use it on jupyter notebook
#this assigment can be found on : /projects/bgmp/dchilin/Bi624/Demultiplex
#copy the files from talapas to this directory using cp
# $ cp /projects/bgmp/shared/2017_sequencing/*gz .

#there should be 4 fastq files
#from the previous demultiplex assignment in Bi622, we determined that: 
    #Read 1:	1294_S1_L008_R1_001.fastq.gz 
	#Index 1	1294_S1_L008_R2_001.fastq.gz 
	#Read 2: 	1294_S1_L008_R4_001.fastq.gz 
	#Index 2: 	1294_S1_L008_R3_001.fastq.gz 
#to verify, less the file and look at the length of the sequences. The short seq (6-10) is the index. 

#create test files 
# $ zcat 1294_S1_L008_R4_001.fastq.gz | head -500 >R4_test.fq
#each file has 500 lines

#download test files into local computer so you can work on jupyter notebook
# $ scp t1:/projects/bgmp/dchilin/Bi624/Demultiplex/*.fq .


###########################################################
#files:
#create variables for each file

#original files
R1 = '/projects/bgmp/dchilin/Bi624/Demultiplex/1294_S1_L008_R1_001.fastq'
I1 = '/projects/bgmp/dchilin/Bi624/Demultiplex/1294_S1_L008_R2_001.fastq'
R2 = '/projects/bgmp/dchilin/Bi624/Demultiplex/1294_S1_L008_R4_001.fastq'
I2 = '/projects/bgmp/dchilin/Bi624/Demultiplex/1294_S1_L008_R3_001.fastq'


#test files contain 1 million lines each
#R1 = '/Users/daisychilin/shell/Bi624/Demultiplex/R1_test.fq'
#I1 = '/Users/daisychilin/shell/Bi624/Demultiplex/R2_test.fq'
#R2 = '/Users/daisychilin/shell/Bi624/Demultiplex/R4_test.fq'
#I2 = '/Users/daisychilin/shell/Bi624/Demultiplex/R3_test.fq'


#############################################################
#FUNCTIONS:

#define functions to convert letters to phred scores
def convert_phred(letter):
    phred_score = ord(letter)-33
    return(phred_score)


#define a function where we can reverse the sequence of an index
def reverse_complement(bases):
    comp_bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return ''.join([comp_bases[base] for base in bases[::-1]])

###############################################################

#create a dictionary for the barcodes
barcodes = {'TAGCGTA':'1','CGATCGAT':'2','GATCAAGG':'3','AACAGCGA':'4', \
            'TAGCCATG':'5', 'CGGTAATC':'6','CTCTGGAT':'7','TACCGGAT':'8', \
            'CTAGCTCA':'9','CACTTCAC':'10','GCTACTCT':'11','ACGATCAG':'12', \
            'TATGGCAC':'13', 'TGTTCCGT':'14','GTCCTAAG':'15', 'TCGACAAG':'16', \
            'TCTTCGAC':'17', 'ATCATGCG':'18', 'ATCGTGGT':'19', 'TCGAGAGT':'20', \
            'TCGGATTC':'21', 'GATCTTGC':'22','AGAGTCCA':'23', 'AGGATAGC':'24'}                                                                                
###############################################################

#counters:
mapped_barcodes = 0
unmapped_barcodes = 0
N_counter_R1 = 0
N_counter_R2 = 0
low_QS_R1 = 0
low_QS_R2 = 0
good_reads = 0


################################################################

#open all 4 files at the same time, use 'r' to read
#open 'trash' files, use 'w' to write into it
with open (R1, 'r') as read_1, \
    open (I1, 'r') as index_1,\
    open (R2, 'r') as read_2, \
    open (I2, 'r') as index_2, \
    open('undetermine_R1.fq', 'w') as un_r1, \
    open('undetermine_R2.fq', 'w') as un_r2:
    while R1 and I2 and R2 and I2:
        #Read1 (R1)
        #remove the new line from each line using .strip()
        header_r1 = read_1.readline().strip('\n')
        if not header_r1:
            break
        seq_r1 = read_1.readline().strip('\n')
        plus_r1 = read_1.readline().strip('\n')
        qs_r1 = read_1.readline().strip('\n')
        #this variable is to combine all the lines together to create our fastq file
        r1_lines = header_r1+'\n' + seq_r1+'\n' + plus_r1+'\n' + qs_r1+'\n'
        
        #Read1 index (I1)
        header_i1 = index_1.readline().strip('\n')
        if not header_i1:
            break
        seq_i1 = index_1.readline().strip('\n')
        plus_i1 = index_1.readline().strip('\n')
        qs_i1 = index_1.readline().strip('\n')
        #print(seq_i1)
        #checking for quality score
        #creat an empty list to input your quality scores in there     
        list_qs_r1 = []
        for charater in qs_i1:    
            score = convert_phred(charater) 
            #append each score to the empty list using .append() 
            list_qs_r1.append(score)
        #sum the list of quality scores
        qs_sum = sum(list_qs_r1)
        #get average by dividing the sum with the total of numbers in line . Use len() 
        avg_r1 = (qs_sum/len(qs_i1))
        #print(avg_r1)
        #set quality score cutoff, we only want reads that meet a certain quality score
        #I chose 30 because a qs of 30 is 99.9% correct. Above 40 is 99.99% . 
        #Both  indexes seem to range from 30-38 which is a hight quality score. 
        #If I cutoff anything above 30, more than 50% of my data will get cut off.
        QS_cutoff = 30 
        if avg_r1 < QS_cutoff:
            #out put this into the 'trash' file
            #keep track of how many reads are going into this trash file
            low_QS_R1 += 1
            un_r1.write(r1_lines)
        else:
            #determining Ns in the index 
            if 'N' in seq_i1:
                #if we see an N present in the index, we are tossing that into the 'trash' file and keeping count of that 
                N_counter_R1 += 1
                un_r1.write(r1_lines)     
                
       
        
        
        #Read2 (R2)
        header_r2 = read_2.readline().strip('\n')
        if not header_r2:
            break
        seq_r2 = read_2.readline().strip('\n')
        plus_r2 = read_2.readline().strip('\n')
        qs_r2 = read_2.readline().strip('\n')
        r2_lines = header_r2+'\n' + seq_r2+'\n' + plus_r2+'\n' + qs_r2+'\n'
        
        #Read2 index (I2)
        header_i2 = index_2.readline().strip('\n')
        if not header_i2:
            break
        seq_i2 = index_2.readline().strip('\n')
        plus_i2 = index_2.readline().strip('\n')
        qs_i2 = index_2.readline().strip('\n')
        #checking for quality score
        list_qs_r2 = []
        for charater in qs_i2:    #for every character(x) in the phred score string
            score_r2 = convert_phred(charater) #defining score as every chacter(x) intro a a converting into a number
            #print(score)
            list_qs_r2.append(score_r2)
        qs_sum_r2 = sum(list_qs_r2)
        #print(qs_sum)
        avg_r2 = (qs_sum_r2/len(qs_i2))
        #print(avg)  
            #create a variable for a quality score cutoff
        QS_cutoff = 33
        if avg_r2 < QS_cutoff:
            low_QS_R2 += 1
            un_r2.write(r2_lines)
        else:
            #determining Ns in the index 
            if 'N' in seq_i2:    
                N_counter_R2 += 1
                #out put the matching lines from R1 to the "trash" file
                un_r2.write(r2_lines)     
            else:
                #take the reverse complement of the index in read 2 
                reverse_comp = reverse_complement(seq_i2)
                #print(reverse_comp)
                if reverse_comp not in barcodes:
                    #if this is not found in the barcode dictionary, toss this read into the 'trash' file 
                    unmapped_barcodes += 1
                    un_r1.write(r1_lines)
                    un_r2.write(r2_lines)
                    #if not(else), just make it equal to one   
                else:  
                    mapped_barcodes += 1
                    if reverse_comp == seq_i1:
                        #if they do match, input this into indvidual files 
                        good_reads += 1
                        #creates a new file using the barcode , but with the .fq extension
                        #use 'a' to append data into the file
                        f_R1 = open(reverse_comp + "_R1.fq", "a") 
                        f_R2 = open(reverse_comp + "_R2.fq", "a") 
                        f_R1.write(r1_lines)
                        f_R2.write(r2_lines)
                     #indexes that did not match (i2 != i1), output them into the 'trash' files
                    else:
                        un_r1.write(r1_lines)
                        un_r2.write(r2_lines)



print('# of reads with unacceptable low quality - read 1:', low_QS_R1)   
print('# of reads with unacceptable low quality - read 2:', low_QS_R2)
print('# of undertermined base calls (N)- read 1:', N_counter_R1)                   
print('# of undertermined base calls (N) - read 2:', N_counter_R2)                   
print('# of unacceptable index pairs:', unmapped_barcodes)
print('# of indexes that match to a barcode:', mapped_barcodes)
print('# of acceptable index pair reads:', good_reads)



#################################################################

#output files:

# AACAGCGA_R1.fq		AACAGCGA_R2.fq
# ACGATCAG_R1.fq		ACGATCAG_R2.fq
# AGAGTCCA_R1.fq 		AGAGTCCA_R2.fq
# AGGATAGC_R1.fq		AGGATAGC_R2.fq
# ATCATGCG_R1.fq		ATCATGCG_R2.fq
# ATCGTGGT_R1.fq		ATCGTGGT_R2.fq
# CACTTCAC_R1.fq		CACTTCAC_R2.fq
# CGATCGAT_R1.fq		CGATCGAT_R2.fq
# CGGTAATC_R1.fq		CGGTAATC_R2.fq
# CTAGCTCA_R1.fq		CTAGCTCA_R2.fq
# CTCTGGAT_R1.fq		CTCTGGAT_R2.fq
# GATCAAGG_R1.fq		GATCAAGG_R2.fq
# GATCTTGC_R1.fq		GATCTTGC_R2.fq
# GCTACTCT_R1.fq		GCTACTCT_R2.fq
# GTCCTAAG_R1.fq		GTCCTAAG_R2.fq
# TACCGGAT_R1.fq		TACCGGAT_R2.fq
# TAGCCATG_R1.fq		TAGCCATG_R2.fq
# TATGGCAC_R1.fq		TATGGCAC_R2.fq
# TCGACAAG_R1.fq		TCGACAAG_R2.fq
# TCGAGAGT_R1.fq		TCGAGAGT_R2.fq
# TCGGATTC_R1.fq		TCGGATTC_R2.fq
# TCTTCGAC_R1.fq		TCTTCGAC_R2.fq
# TGTTCCGT_R1.fq		TGTTCCGT_R2.fq
# undetermine_R1.fq	undetermine_R2.fq




##################################################################        
        
# test files
#R1 = '/Users/daisychilin/shell/Bi624/Demultiplex/test_files/test_R1.fq'
#I1 = '/Users/daisychilin/shell/Bi624/Demultiplex/test_files/test_I1.fq'
#R2 = '/Users/daisychilin/shell/Bi624/Demultiplex/test_files/test_R2.fq'
#I2 = '/Users/daisychilin/shell/Bi624/Demultiplex/test_files/test_I2.fq'

#barcodes = {'AAAAA':'1','TTTTT': '2', 'CCCCC':'3', 'GGGGG':'4' }                                                                                    

#note: 
#each file has 3 reads each. We should expect: . 
    # 1 read with bad QS
    # 1 read with undertermined baes call (N)
    # 1 good read
    
#Results:
# reads found with low quality score -R1: 1
# reads found with low quality score -R2: 1
# undertermined base calls (N)- R1: 1
# undertermined base calls (N) - R2: 1
# indexes that match to a barcode: 1
# indexes that did not match to a barcode: 0
# demultiplexed reads: 1 

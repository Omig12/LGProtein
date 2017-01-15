##########################################################################################
#                                                                                        #
# Copyright 2016 Megaprobe-Lab                                                           #
#                                                                                        #
# This is software created by the megaprobe lab under the GPL3 license.                  #
#                                                                                        #
# This toy program parses genetic sequences and produces an equivalent list of the       #
# amino acids contained within the sequence. To run the program, utilize the following   #
# command: "Python2.7 n2aa.py [FILE | Sequence]". The FILE parameter should              #
# contain a file with at least one genetic sequence (DNA or RNA).                        #
#                                                                                        #
##########################################################################################

# Load Modules
import os
import sys

############################################
#                                          #
# Read the sequence and translate it       #
#                                          #
############################################

def parse(seq):

    # Amino Acids Dicitonary 
    DIC = {'ALA': ['GCT','GCC','GCA', 'GCG'], 'GLY': ['GGT', 'GGC', 'GGA', 'GGG'], 'ILE':['ATT', 'ATC', 'ATA'], 'LEU':['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'], 'PRO': ['CCT', 'CCC', 'CAA', 'CCG'], 'VAL' : ['GTT', 'GTC', 'GTA', 'GTG'], 'PHE': ['TTT','TTC'], 'TRP': ['TGG'], 'TYR': ['TAT', 'TAC'], 'ASP': ['GAT', 'GAC'], 'GLU': ['GAA', 'GAG'],  'ARG': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'HIS': ['CAT', 'CAC'], 'LYS': ['AAA', 'AAG'], 'SER': ['TCT','TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'THR': ['ACT', 'ACC', 'ACA', 'ACG'], 'CYS': ['TGT', 'TGC'], 'MET': ['ATG'], 'ASN' : ['AAT', 'AAC'], 'GLN': ['CAA', 'CAG']}
    
    single = {'ALA': 'A', 'GLY': 'G', 'ILE':'I', 'LEU':'L', 'PRO': 'P', 'VAL' : 'V', 'PHE': 'F', 'TRP': 'W', 'TYR': 'Y', 'ASP': 'D', 'GLU': 'E',  'ARG': 'R', 'HIS': 'H', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'CYS': 'C', 'MET': 'M', 'ASN' : 'N', 'GLN': 'Q'}
    
    # Convert RNA to DNA 
    seq = seq.upper().strip().replace('U', 'T')
    
    # Validate Input
    for i in seq:
        try:
            if i not in 'ATCG':
                raise ValueError ('Illegal Character')
        except ValueError:    
            print ('Invalid character "'+ i +'" in input.')
            raise 

    # Unfairly trim sequence to required length
    seq =  seq[0:len(seq)-(len(seq) % 3)]

    # Array to store translated sequence 
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]


    # Array to hold amino acids 
    AA = []
    
    # Translate
    for i in codons: 
        for key, value in DIC.iteritems():
            if i in value:
                AA.append(key)

############################################
#                                          #
# Start writing output of translation      #
#                                          #
############################################	
    
    if ('./amino-acids.txt'  exists ):
    with open('./amino-acids.txt', 'w+') as output:
        output.write('The sequence: \n\t'+ seq +'\n\nTraslates to: \n\t') 
        for i in AA:   
            output.write(i + '-')
        
        output.write('\n\nSingle letter format: \n\t')
        for i in AA:
        	output.write(single[i]+ '-') 

###################
#                 #
#  Read input     #
#                 #
###################
# def get_input():
# Initiate file processing
try:  
    arg = sys.argv[1]

    try: 
        arr = ""
        with open(arg, 'r') as f:
            for line in f:
                arr += line.strip()
        arg = arr

    except:
        pass

except:
    arg = raw_input("Please input sequence: ")

finally: 
    parse(arg)
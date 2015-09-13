# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:21:46 2015

@author: tyler
"""

from Bio import SeqIO
import os

def main():
    p_path = "./plasmids/"
    for dirName, subdirList, fileList in os.walk(p_path):
        for f in fileList:
            my_seq = SeqIO.read(p_path+f, "genbank")
            annotation_check(my_seq)
        
    return
    
    # We need to check all relevant annotations
    # prefix, phage origins, promoter, rbs, cds, suffix
def annotation_check(my_seq):
    print "Checking sequence", my_seq.id + ":"
    path = "./sequences/" # location of correct part sequences
    promoter_check = False # will be set to true if a promoter annotation exists
    rbs_check = False
    cds_check = False
    seq_dict = ({'tactag':'tactag', 'tactagag':'tactagag', 'dna: biobrick prefix':'attcgcggccgcttctagag',
                 'dna: biobrick suffix':'tactagtagcggccgctgcag'})
    error_list = list()
    my_seq_nts = my_seq.seq
    for feat in my_seq.features:
        feat_type = feat.type # rbs, promoter, etc
        
        if ("RBS" in feat_type): rbs_check = True
        if ("promoter" in feat_type): promoter_check = True
        if ("CDS" in feat_type): cds_check = True
        seq1_indices = feat.location # location on plasmid            
        start = int(seq1_indices.start) # start index of feature
        end = int(seq1_indices.end) # end index of feature
        feat_seq = my_seq_nts[start:end].lower()
        for key in feat.qualifiers:              

            feat_key = (feat.qualifiers[key][0]).lower() # go through every key in feature  
            if (os.access(path+feat_key+".gb", os.F_OK)):
                seq2 = next(SeqIO.parse(path+feat_key+".gb", "genbank"))
                if (feat_key == seq2.id.lower()): # this feature is one of our biobricks
                    if (str(feat_seq) == str(seq2.seq.lower())):
                        continue
                    else:
                        error_list.append(seq2.id)
            elif (feat_key in seq_dict): # check other features
                seq2 = seq_dict[feat_key]
                if (str(feat_seq) == str(seq2)):
                    continue
                else:
                    error_list.append(feat_key)
    if not(rbs_check): error_list.append("Missing RBS")
    if not(promoter_check): error_list.append("Missing promoter")
    if not(cds_check): error_list.append("Missing CDS")
    
    if (len(error_list) == 0):
        print "All annotations valid.\n"
    else:
        print "Invalid annotations:"
        for error in error_list:
            print error
    return
    
main()
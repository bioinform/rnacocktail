#!/usr/bin/python

############################################################################
#This script is modified from the original code by Kin Fai Au 
#Obtained from https://github.com/jason-weirather/Au-public/blob/master/gold/gpd2gtf.py
#Available unde Apache License Version 2.0
############################################################################

import sys
import math

### generate_transcript_list
############################
def generate_transcript_list(gpd_file, transcript_list):
    
    for line in gpd_file:
        
        if (line[0] == '#'):
            continue
        
        fields = line.split()
        num_exons = int(fields[8])

        start_pos_list = fields[9].split(',')
        end_pos_list = fields[10].split(',')
        
        exon_pos = [0] * num_exons
        for i in range(num_exons):
           exon_pos[i] = [start_pos_list[i], end_pos_list[i]]
        
        transcript_list.append([fields[0], fields[1], fields[2], fields[3], exon_pos])
        

### generate_FPKM_dict
#######################
def generate_FPKM_dict(FPKM_file, FPKM_dict):
    
    for line in FPKM_file:
        fields = line.split()
        FPKM_dict[fields[0]] = fields[1]
        
    

### generate_gpd_format
#######################
def generate_gtf_format(gtf_file, transcript_list, FPKM_dict, source):

    
    for line in transcript_list:
        exon_pos = line[4]
        # transcript line
        
        # chr name 
        gtf_file.write(line[2] + '\t' + source + '\t' + "transcript" + '\t')
        # start-end pos, score
        gtf_file.write("%s"%(int(exon_pos[0][0])+1) + '\t' + exon_pos[-1][1] + '\t' + '*' + '\t')
        # Direction
        gtf_file.write(line[3] + '\t' + '.' + '\t')
        
        if (FPKM_dict.has_key( line[1]) ):
            FPKM = FPKM_dict[line[1]]
        else:
            FPKM = '*'
        attribute_1 = 'gene_id "' + line[0] + '"; transcript_id "' + line[1] + '"; '
        attribute_2 = ('FPKM "' + FPKM + '"; frac "' + '*' + '"; conf_lo "' + '*' + '"; ' +
                       'conf_hi "' + '*' + '"; cov "' + '*' + '";\n')
        
        gtf_file.write(attribute_1)
        gtf_file.write(attribute_2)
        
        num_exons = len(exon_pos)
        for i in range(num_exons):
            # chr name 
            gtf_file.write(line[2] + '\t' + source + '\t' + "exon" + '\t')
            # start-end pos, score
            gtf_file.write("%s"%(int(exon_pos[i][0])+1) + '\t' + exon_pos[i][1] + '\t' + '*' + '\t')
            # Direction
            gtf_file.write(line[3] + '\t' + '.' + '\t')
            gtf_file.write(attribute_1)
            gtf_file.write('exon_number "' +  str(i+1) + '"; ')
            gtf_file.write(attribute_2)



### Main
########
def main():
    gpd_file = open(sys.argv[1], 'r')
    FPKM_file = open(sys.argv[2], 'r')
    gtf_file = open(sys.argv[3], 'w')
    source = sys.argv[4]

    transcript_list = []
    FPKM_dict = dict()
    generate_transcript_list(gpd_file, transcript_list)
    generate_FPKM_dict(FPKM_file, FPKM_dict)
    generate_gtf_format(gtf_file, transcript_list, FPKM_dict, source)
    
    gpd_file.close()
    gtf_file.close()


if __name__ == '__main__':
    main()

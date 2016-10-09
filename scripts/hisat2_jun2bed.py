#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2:
    HISATJun_filename = sys.argv[1]
    bed_filename = sys.argv[2]
else:
    print("usage: python hisat2_jun2bed.py HISAT2_splicesites.txt junction.bed")
    sys.exit(1)

jun_s = set()
junction=open(HISATJun_filename,'r')
output_file = open(bed_filename,'w')
for line in junction:
    if line[0:5]=='track':
        continue
    else:
        line_list=line.strip().split("\t")
        leftpos=str(int(line_list[1]))
        rightpos=str(int(line_list[2]))
        locus = "_".join([line_list[0],leftpos,rightpos,line_list[3]])
        jun_s.add(locus)

output_file.write("track name=junctions description=\"HISAT2 junctions\"\n")
i=0
for locus in jun_s:
    output_ls = []
    locus_ls = locus.split("_")
    chr_name = locus_ls[0]
    int_start =  int(locus_ls[1])-51
    if int_start<=0:
        start = "1"
        width_start = str(49+int_start)
    else:
        start = str(int_start)
        width_start = "50"
    end = str( int(locus_ls[2]) + 50 )
    distance = str( int(locus_ls[2])  - int(locus_ls[1])+51 )

    sign = locus_ls[3]

    name = "HISAT" + str(i)

    i += 1
    output_ls = [chr_name,start,end,name,"50",sign,start,end,"0,0,0","2",width_start+",50","0,"+distance]
    output_file.write( '\t'.join(output_ls) + "\n" )
junction.close()
output_file.close()


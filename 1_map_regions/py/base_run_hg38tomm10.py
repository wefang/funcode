#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 02:57:12 2022

@author: celeste
"""

import base_module as bm
from collections import defaultdict

inputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/hg38_mm10/ENCFF503GCK.tsv'
outputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/new_alignment_result/hg38tomm10_base_alignment.txt'
o = open(outputfile, 'w')
o.write('identifier\thg38_summit\tmm10_summit\thg38_extend_region\tmm10_extend_region\taxtnet_id\taxtnet_direction\tgap_status\n')
error = []

# coordiantes transformation if alignment on negative strand
rev_spec = 'mm10'
#load chromosome information
dict_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/new_align_script_github/' + rev_spec + '.chrom.sizes.txt'
rev_chrom_dict = bm.chrom_size(dict_fp)
# axt net file suffix
axt_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/hg38_mm10/axtnet/'
axt_suffix = '.hg38.mm10.net.axt'
# extent region length
setlen = 200

chrom_0 = 'chr0'
with open(inputfile,'r') as i:
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[3]
            chrom = linelist[0]
            base = int(linelist[6]) + 1
            refile = axt_fp + chrom + axt_suffix
            
            if chrom != chrom_0:
                block_index = 0
                block_dict =  defaultdict()           
                for each in bm.readblock(refile):
                    block_dict[block_index] =  each[:]
                    block_index = block_index + 1
                    
            posi = bm.binarySearch(0, len(block_dict)-1, base, block_dict)
            if posi != 'no match':
                out = bm.outputline(posi, base, setlen, block_dict, rev_chrom_dict)
                final_line = index + '\t' + out
                o.write(final_line)
            chrom_0 = chrom
        except Exception as e:
            print(e)
            print(line)
            error.append(index)
o.close()
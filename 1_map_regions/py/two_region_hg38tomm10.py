#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 16:19:39 2022

@author: celeste
"""

import region_module as rm
from collections import defaultdict

inputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/hg38_mm10/hg38_mm10_region.txt'
outputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/hg38_mm10_phast_phylop/hg38tomm10_keep_hg_percentage_0.txt'
o = open(outputfile,'w')
error = []

# coordiantes transformation if alignment on negative strand
rev_spec = 'mm10'
#load chromosome information of aligned species
dict_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/new_align_script_github/' + rev_spec + '.chrom.sizes.txt'
rev_chrom_dict = rm.chrom_size(dict_fp)
# axt net file suffix
axt_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/hg38_mm10/axtnet/'
axt_suffix = '.hg38.mm10.net.axt'

chrom_0 = 'chr0'
with open(inputfile,'r') as i:
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[1]
            ref_region = linelist[2]
            cmp_region = linelist[3]
            chrom = ref_region.split(':')[0]
            start = int(ref_region.split(':')[1].split('-')[0])
            end = int(ref_region.split(':')[1].split('-')[1]) 
            mylist = [int(start),int(end)]
            cmplist = [int(cmp_region.split(':')[1].split('-')[0]), int(cmp_region.split(':')[1].split('-')[1])]
            refile = axt_fp + chrom + axt_suffix 
            
            if chrom != chrom_0:
                block_index = 0
                block_dict =  defaultdict()           
                for each in rm.readblock(refile):
                    block_dict[block_index] =  each[:]
                    block_index = block_index + 1
            
            con,posi = rm.binarySearch(0, len(block_dict)-1, mylist, block_dict)
            if con == 'in':
                out = rm.outputline_ref_percentage([posi], mylist, block_dict, rev_chrom_dict)
                if out != 'no match\n':
                    final_line = index + '\t' + out
                    o.write(final_line)
                    
            if con in ['ovlow','ovhigh','out']:
                posilist = rm.contSearch(con, posi, mylist, block_dict)
                out = rm.outputline_ref_percentage(posilist, mylist, block_dict, rev_chrom_dict)
                if out != 'no match\n':
                    final_line = ''
                    for eachline in out.split('\n'):
                        if eachline != '':
                            final_line = final_line + index + '\t' + eachline + '\n'
                    o.write(final_line)
            chrom_0 = chrom
            
        except Exception as e:
            print(e)
            print(line)
            error.append(index)
o.close()
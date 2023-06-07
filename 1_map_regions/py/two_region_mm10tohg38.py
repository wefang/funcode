#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 23:26:56 2022
@author: celeste
"""

import region_module as rm
from collections import defaultdict

# inputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/testfile.txt'
# outputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/testout.txt'

inputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/mm10tohg38_base_alignment.txt'
outputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38_phast_phylop/mm10tohg38_keep_hg_percentage_0.txt'
o = open(outputfile,'w')
error = []

# coordiantes transformation if alignment on negative strand
rev_spec = 'hg38'
#load chromosome information of aligned species
dict_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/new_align_script_github/' + rev_spec + '.chrom.sizes.txt'
rev_chrom_dict = rm.chrom_size(dict_fp)
# axt net file suffix
axt_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/axtnet/'
axt_suffix = '.mm10.hg38.net.axt'

chrom_0 = 'chr0'
with open(inputfile,'r') as i:
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[0]
            ref_region = linelist[3]
            cmp_region = linelist[4]
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
                out = rm.outputline_cmp_percentage([posi], cmplist, block_dict, rev_chrom_dict)
                if out != 'no match\n':
                    final_line = index + '\t' + out
                    o.write(final_line)
                    
            if con in ['ovlow','ovhigh','out']:
                posilist = rm.contSearch(con, posi, mylist, block_dict)
                out = rm.outputline_cmp_percentage(posilist, cmplist, block_dict, rev_chrom_dict)
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

 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 23:42:25 2022

@author: celeste
"""
# generator of each alignment block in axtnet file
def readblock(fp):
    lines = []
    with open(fp) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line != '\n':
                lines.append(line.strip())
            if line == '\n':
                yield lines
                lines.clear()
        if lines!= []:
            yield line

def chrom_size(fp):
    mydict = {}
    with open(fp, 'r') as f:
        for line in f:
            (chrom, size) = line.strip().split()
            mydict[chrom] = int(size)
    return(mydict)

def get_ture_index(seq):
    true_index = []
    current = -1
    for i in seq:
        if i != '-':
            current += 1
        true_index.append(current)
    return true_index


def seqdict(ref_seq,cmp_seq):
    seq1 = list(ref_seq.strip().upper())
    seq2 = list(cmp_seq.strip().upper())
    index1 = get_ture_index(seq1)
    index2 = get_ture_index(seq2)
    
    seq_dict = {}
    for i in range(len(seq1)):
        seq_dict[i] = {'ref_trueindex': index1[i], 'cmp_trueindex': index2[i], 'ref_base': seq1[i], 'cmp_base': seq2[i]}  
    return seq_dict


def searchclosest(key, seqdict):
    start = key
    end = key
    while start >= 0:
        if seqdict[start]['cmp_base'] != '-':
            mystart = start
            break
        start = start - 1
    while end <= len(seqdict)-1:
        if seqdict[end]['cmp_base'] != '-':
            myend = end
            break
        end = end + 1
    inter1 = key - mystart
    inter2 = myend - key
    if inter1 <= inter2:
        return 'before', mystart
    if inter1 > inter2:
        return 'after', myend
       
def givcon(midlist,base):
    midstart = int(midlist[0])
    midend = int(midlist[1])
    mybase = int(base)
    if mybase < midstart:
        return 'low'
    if mybase > midend:
        return 'high'
    if mybase >= midstart and mybase <= midend:
        return 'in'

def binarySearch(low, high, base, blockdict):
    if high >= low:
        mid = int(low + (high - low)/2)
        mid_list = blockdict[mid][0].split()[2:4]
        if givcon(mid_list,base) == 'low':
            return(binarySearch(low, mid-1, base, blockdict))
        if givcon(mid_list,base) == 'high':
            return(binarySearch(mid+1, high, base, blockdict))
        if givcon(mid_list,base) == 'in':
            return mid
    else:
        return 'no match'

def extend_region(base, setlen, chromlen):
    thestart = base - int(setlen/2) + 1
    thend = base + int(setlen/2)
    if thestart <= 0:
        thestart = 1
    if thend > chromlen:
        thend = chromlen
    return str(thestart), str(thend)

def outputline(posi, base, setlen, blockdict, ref_chrom_dict, cmp_chrom_dict):
    theline = blockdict[posi][0].split()
    align_id, ref_chrom, ref_start, cmp_chrom, cmp_start, align_dire = theline[0],theline[1],int(theline[2]),theline[4],int(theline[5]),theline[7]
    seq_dict = seqdict(blockdict[posi][1],blockdict[posi][2])
    ref_chrom_len = ref_chrom_dict[ref_chrom]
    cmp_chrom_len = cmp_chrom_dict[cmp_chrom]
    
    for key,value in seq_dict.items():
        if ref_start + value['ref_trueindex'] == base:
            out_cmp = cmp_start + value['cmp_trueindex']
            if value['cmp_base'] != '-':
                status = 'non-gap'
            else:
                status = 'gap'
                closest, thekey = searchclosest(key, seq_dict)
                if closest == 'after':
                    out_cmp = cmp_start + value['cmp_trueindex'] + 1
            break
                    
    if align_dire == '-':
        out_cmp = cmp_chrom_len - out_cmp + 1
    
    my_base = ref_chrom + ':' + str(base) 
    out_posi = cmp_chrom + ':' + str(out_cmp)
    my_region= ref_chrom + ':' + extend_region(base, setlen, ref_chrom_len)[0] + '-' + extend_region(base, setlen, ref_chrom_len)[1]
    out_region = cmp_chrom + ':' + extend_region(out_cmp, setlen, cmp_chrom_len)[0] + '-' + extend_region(out_cmp, setlen, cmp_chrom_len)[1]
    output_line = my_base + '\t' + out_posi + '\t' + my_region + '\t' +  out_region + '\t'  + str(align_id) + '\t'+ align_dire +'\t'+ status + '\n'

    return output_line

   
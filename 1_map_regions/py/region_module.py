#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 21:55:25 2022

@author: celeste
"""

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

def calcoor(chrom_len, start, end):
    tran_start = chrom_len - end + 1
    tran_end = chrom_len - start + 1
    return tran_start, tran_end

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

def get_same_num(start_key, end_key, seqdict):
    same_num = 0
    my_ref_seq = []
    my_cmp_seq = []
    for key, value in seqdict.items():
        if key >= start_key and key <= end_key:
            my_ref_seq.append(value['ref_base'])
            my_cmp_seq.append(value['cmp_base'])
            if value['ref_base'] == value['cmp_base'] and value['ref_base'] != '-':
                same_num = same_num + 1           
    my_ref = ''.join(my_ref_seq)
    my_cmp = ''.join(my_cmp_seq)
    return my_ref, my_cmp, same_num
       
def givcon(midlist,mylist):
    midstart = int(midlist[0])
    midend = int(midlist[1])
    mystart = int(mylist[0])
    myend = int(mylist[1])
    if myend <= midstart:
        return 'low'
    if mystart >= midend:
        return 'high'
    if myend > midstart and myend < midend and mystart < midstart:
        return 'ovlow'
    if mystart > midstart and mystart < midend and myend > midend:
        return 'ovhigh'
    if mystart >= midstart and myend <= midend:
        return 'in'
    if mystart <= midstart and myend >= midend:
        return 'out'
    

def binarySearch(low, high, region, blockdict):
    if high >= low:
        mid = int(low + (high - low)/2)
        mid_list = blockdict[mid][0].split()[2:4]
        if givcon(mid_list,region) == 'low':
            return(binarySearch(low, mid-1, region, blockdict))
        if givcon(mid_list,region) == 'high':
            return(binarySearch(mid+1, high, region, blockdict))
        if givcon(mid_list,region) in ['ovlow','ovhigh','in','out']:
            return givcon(mid_list,region),mid
    else:
        return 'no match', -1



def contSearch(con, midblock, region, blockdict):
    result = [midblock]
    if con == 'ovlow':
        i =  midblock - 1
        while i >= 0 and i < len(blockdict):
            linelist = blockdict[i][0].split()[2:4]
            if givcon(linelist,region) in ['ovlow','ovhigh','in','out']:
                result.append(i)
            else:
                break
            i = i - 1
    if con == 'ovhigh':
        j = midblock + 1
        while j >= 0 and j < len(blockdict):
            linelist = blockdict[j][0].split()[2:4]
            if givcon(linelist,region) in ['ovlow','ovhigh','in','out']:
                result.append(j)
            else:
                break
            j = j + 1
    if con == 'out':
        a = midblock - 1
        b = midblock + 1
        while a >= 0 and a < len(blockdict) :
            lista = blockdict[a][0].split()[2:4]
            if givcon(lista,region) in ['ovlow','ovhigh','in','out']:
                result.append(a)
            else:
                break
            a = a - 1
        while b >= 0 and b < len(blockdict):
            listb = blockdict[b][0].split()[2:4]
            if givcon(listb,region) in ['ovlow','ovhigh','in','out']:
                result.append(b)
            else:
                break
            b = b + 1
    return result



def outputline_ref_percentage(posilist, region, blockdict, rev_chrom_dict):
    id_list = []
    myposi_list = []
    outposi_list = []
    number_list = []
    direct_list = []
    
    for i in posilist:
        theline = blockdict[i][0].split()
        align_id, ref_chrom, ref_start, cmp_chrom, cmp_start, align_dire = theline[0],theline[1],int(theline[2]), theline[4], int(theline[5]), theline[7]
        seq_dict = seqdict(blockdict[i][1],blockdict[i][2])
        my_start, my_end = int(region[0]), int(region[1])
        rev_chrom_len = rev_chrom_dict[cmp_chrom]
        
        key_start, key_end = [], [len(seq_dict) - 1]
        for key,value in seq_dict.items():
            if ref_start + value['ref_trueindex'] == my_start:
                key_start.append(key)
            if ref_start + value['ref_trueindex'] == my_end:
                key_end.append(key)
        if key_start == []:
            key_start = [0]
        key_start, key_end = min(key_start), min(key_end)
                   
        out_ref_start = ref_start + seq_dict[key_start]['ref_trueindex']
        out_ref_end = ref_start + seq_dict[key_end]['ref_trueindex']
        out_cmp_start = cmp_start + seq_dict[key_start]['cmp_trueindex']
        out_cmp_end = cmp_start + seq_dict[key_end]['cmp_trueindex']
        
        if align_dire == '-':
            out_cmp_start, out_cmp_end = calcoor(rev_chrom_len, out_cmp_start, out_cmp_end)
        
        my_posi = ref_chrom + ':' + str(out_ref_start) + '-' + str(out_ref_end)
        out_posi = cmp_chrom + ':' + str(out_cmp_start) + '-' + str(out_cmp_end)
        out_refseq, out_cmpseq, out_number = get_same_num(key_start, key_end, seq_dict)
        
        id_list.append(align_id)
        myposi_list.append(my_posi)
        outposi_list.append(out_posi)
        number_list.append(out_number)
        direct_list.append(align_dire)
            
    if id_list == []:
        output_line = 'no_match\n'
    else:
        seq_length = int(region[1]) - int(region[0]) + 1
        output_line = ''
        for k in range(0,len(id_list)):
            output_line = output_line + myposi_list[k] + '\t' + outposi_list[k]+'\t'+ str(id_list[k])+'\t'+\
                            str(direct_list[k])+'\t'+ str(number_list[k]) + '\t' + format(sum(number_list)/seq_length*100,'0.3f') + '\n'
    return output_line


def outputline_cmp_percentage(posilist, cmp_region, blockdict, rev_chrom_dict):
    id_list = []
    myposi_list = []
    outposi_list = []
    number_list = []
    direct_list = []
    
    for i in posilist:
        theline = blockdict[i][0].split()
        align_id, ref_chrom, ref_start, cmp_chrom, cmp_start, cmp_end, align_dire = theline[0],theline[1],int(theline[2]),theline[4],int(theline[5]), int(theline[6]),theline[7]
        seq_dict = seqdict(blockdict[i][1],blockdict[i][2])
        my_start, my_end = cmp_region[0], cmp_region[1]
        rev_chrom_len = rev_chrom_dict[cmp_chrom]
        if align_dire == '-':
            my_start, my_end = calcoor(rev_chrom_len, my_start, my_end)
        if my_start > cmp_end or my_end < cmp_start:
            continue
        
        key_start, key_end = [], [len(seq_dict) - 1]
        for key,value in seq_dict.items():
            if cmp_start + value['cmp_trueindex'] == my_start:
                key_start.append(key)
            if cmp_start + value['cmp_trueindex'] == my_end:
                key_end.append(key)
        if key_start == []:
            key_start = [0]
        key_start, key_end = min(key_start), min(key_end)
        
        
        if key_start >= key_end:
            continue
        else:
            out_ref_start = ref_start + seq_dict[key_start]['ref_trueindex']
            out_ref_end = ref_start + seq_dict[key_end]['ref_trueindex']
            out_cmp_start = cmp_start + seq_dict[key_start]['cmp_trueindex']
            out_cmp_end = cmp_start + seq_dict[key_end]['cmp_trueindex']
                    
            if align_dire == '-':
                out_cmp_start, out_cmp_end = calcoor(rev_chrom_len, out_cmp_start, out_cmp_end)

            
            my_posi = ref_chrom + ':' + str(out_ref_start) + '-' + str(out_ref_end)
            out_posi = cmp_chrom + ':' + str(out_cmp_start) + '-' + str(out_cmp_end)
            out_refseq, out_cmpseq, out_number = get_same_num(key_start, key_end, seq_dict)
    
                
            id_list.append(align_id)
            myposi_list.append(my_posi)
            outposi_list.append(out_posi)
            number_list.append(out_number)
            direct_list.append(align_dire)
            
    if id_list == []:
        output_line = 'no_match\n'
    else:
        seq_length = int(cmp_region[1]) - int(cmp_region[0]) + 1
        output_line = ''
        for k in range(0,len(id_list)):
            output_line = output_line + myposi_list[k] + '\t' + outposi_list[k]+'\t'+ str(id_list[k])+'\t'+\
                            str(direct_list[k])+'\t'+ str(number_list[k]) + '\t' + format(sum(number_list)/seq_length*100,'0.3f') + '\n'
    return output_line


# def outputline_2region(fp, posilist, ref_region, cmp_region, blockdict, rev_spec):
#     id_list = []
#     myposi_list = []
#     outposi_list = []
#     percen_list = []
#     direct_list = []
    
#     for i in posilist:
#         theline = blockdict[i][0].split()
#         align_id, ref_chrom, ref_start, cmp_chrom, cmp_start, align_dire = theline[0],theline[1],int(theline[2]),theline[4],int(theline[5]),theline[7]
#         seq_dict = seqdict(blockdict[i][1],blockdict[i][2])
#         my_start1, my_end1 = ref_region[0], ref_region[1]
#         my_start2, my_end2 = cmp_region[0], cmp_region[1]
#         if align_dire == '-':
#             my_start2, my_end2 = calcoor(rev_spec, cmp_chrom, cmp_region[0], cmp_region[1])
            
#         key_start1, key_end1, key_start2, key_end2 = 0, len(seq_dict) - 1, 0, len(seq_dict) - 1
#         for key,value in seq_dict.items():
#             if ref_start + value[0] == my_start1:
#                 key_start1 = key
#             if ref_start + value[0] == my_end1:
#                 key_end1 = key
#             if cmp_start + value[1] == my_start2:
#                 key_start2 = key
#             if cmp_start + value[1] == my_end2:
#                 key_end2 = key
        
#         key_start = max([key_start1, key_start2])
#         key_end = min([key_end1, key_end2])
        
#         if key_start >= key_end:
#             continue
#         else:
#             out_ref_start = ref_start + seq_dict[key_start][0]
#             out_ref_end = ref_start + seq_dict[key_end][0]
#             out_cmp_start = cmp_start + seq_dict[key_start][1]
#             out_cmp_end = cmp_start + seq_dict[key_end][1]
                    
#             if align_dire == '-':
#                 out_cmp_start, out_cmp_end = calcoor(rev_spec, cmp_chrom, out_cmp_start, out_cmp_end)
#                 align_dire = '+'
            
#             my_posi = ref_chrom + ':' + str(out_ref_start) + '-' + str(out_ref_end)
#             out_posi = cmp_chrom + ':' + str(out_cmp_start) + '-' + str(out_cmp_end)
#             out_refseq, out_cmpseq, out_percent = calpercent(key_start, key_end, seq_dict)
    
                
#             id_list.append(align_id)
#             myposi_list.append(my_posi)
#             outposi_list.append(out_posi)
#             percen_list.append(out_percent)
#             direct_list.append(align_dire)
            
#     if id_list == []:
#         output_line = 'no_match\n'
#     else:
#         output_line = ''
#         for k in range(0,len(id_list)):
#             output_line = output_line + myposi_list[k] + '\t' + outposi_list[k]+'\t'+ str(id_list[k])+'\t'+\
#                             str(direct_list[k])+'\t'+str(percen_list[k]) + '\n'
#     return output_line

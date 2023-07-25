#coding=utf-8

'''
This script generates splicing amount of given circRNAs in samples based on junction provided by annotation file (gtf).
Annotation of version gencode.v38lift37 is recommended.
GTF files of circRNA information generated with CIRIquant and SAM (or BAM) files  generated with BWA are needed.

Usage:
python ./script_splicingamount.py /path/to/inputfile /output/dir /path/to/gtffile
Input file contains three columns seperated by tab. Header is needed. 
The first column contains the sample index, the second column contains the path to CIRIquant derived gtf file and the third column contains the path to SAM file.
'''

import os
from subprocess import *
import re
import sys
import pickle

def gtf_junction_extract(fn):
    '''
    extract junction information from genecode gtf
    return junction_dict
    '''
    CHR=0
    GENE_PART=2
    EXON_START=3
    EXON_END=4
    INFO=8
    
    # initialize junction_dict
    chrs = ['chr' + str(i) for i in range(1, 23)]
    chrs.extend(['chrX', 'chrY', 'chrM'])
    junction_dict = {}
    for c in chrs:
        junction_dict[c] = {'up':{}, 'down':{}}
    
    # read gtf
    with open(fn, 'r') as f:
        for line in f:
            if len(line) == 0:
                continue
            if line.startswith('#'):
                continue
            content = line.strip().split('\t')
            gene_part = content[GENE_PART]
            if gene_part != 'exon':
                continue
            chrom = content[CHR]
            exon_start = content[EXON_START]
            for pos in range(int(exon_start)-6, int(exon_start)+7):
                if pos not in junction_dict[chrom]['down']:
                    junction_dict[chrom]['down'][pos] = {'gene_id':[]}
            exon_end = content[EXON_END]
            for pos in range(int(exon_end)-6, int(exon_end)+7):
                if pos not in junction_dict[chrom]['up']:
                    junction_dict[chrom]['up'][pos] = {'gene_id':[]}
            info = content[INFO]
            gene_id = [x for x in info.split(';') if 'gene_id' in x][0]
            gene_id = gene_id.strip().split(' ')[1].strip('"')
            for pos in range(int(exon_start)-6, int(exon_start)+7):
                if gene_id not in junction_dict[chrom]['down'][pos]['gene_id']:
                    junction_dict[chrom]['down'][pos]['gene_id'].append(gene_id)
            for pos in range(int(exon_end)-6, int(exon_end)+7):
                if gene_id not in junction_dict[chrom]['up'][pos]['gene_id']:
                    junction_dict[chrom]['up'][pos]['gene_id'].append(gene_id)
    return junction_dict

def read_ciri_result(fn):
    '''
    extract exon circ id and correspond gene name and gene id
    return ciri_dict
    '''
    # initiate ciri_dict
    ciri_gene_dict = {}
    ciri_circ_dict = {}
    
    # read ciri_quant result    
    with open(fn, 'r') as f:
        for line in f:
            if len(line) == 0:
                continue
            if line.startswith('#'):
                continue
            content = line.strip().split('\t')
            info = content[-1].strip().split('; ')
            circ_type = [x for x in info if 'circ_type' in x][0]
            circ_type = circ_type.split(' ')[-1].strip('"')
            if circ_type != 'exon':
                continue
            circ_id = [x for x in info if 'circ_id' in x][0]
            circ_id = circ_id.split(' ')[-1].strip('"')
            gene_ids = [x for x in info if 'gene_id' in x][0]
            gene_ids = sorted(gene_ids.split(' ')[-1].strip('"').split(','))
            for gene_id in gene_ids:
                if gene_id not in ciri_gene_dict:
                    ciri_gene_dict[gene_id] = 0 # geneid: splicing amount
            ciri_circ_dict[circ_id] = gene_ids # circid: geneid
    return ciri_gene_dict, ciri_circ_dict

def extract_juncreads_from_sam(sam_fn, ciri_result, junction_dict, outputdir):
        
    # extract juncreads from sam
    MAPQ_INDEX = 4
    CIGAR_INDEX = 5
    CHROM_INDEX = 2
    POS_INDEX = 3
    NAME_INDEX = 0
    FLAG_INDEX = 1
    first_align = False
    align_pattern = None
    align1_gene_id = []
    align2_gene_id = []
    linear_align1 = False
    circ_align1 = False
    linear_align2 = False
    circ_align2 = False
    
    with open(sam_fn, 'r') as f:
        chrs = ['chr' + str(i) for i in range(1, 23)]
        chrs.extend(['chrX', 'chrY', 'chrM'])
        for line in f:
            if len(line) == 0:
                continue
            if line.startswith('@'):
                continue
            content = line.strip().split('\t')
            # MAPQ screening
            MAPQ = content[MAPQ_INDEX]
            if int(MAPQ) <= 5:
                first_align = False
                continue
            # CIGAR screening
            CIGAR = content[CIGAR_INDEX]
            chrom = content[CHROM_INDEX]
            if chrom not in chrs:
                continue
            if first_align:
                align2_gene_id = []
                linear_align2 = False
                circ_align2 = False
                first_align = False
                if align_pattern == 'MS':
                    if re.match(r'^\d+[SH]\d+M$', CIGAR):
                        align2_name = content[NAME_INDEX]
                        align2_chr = content[CHROM_INDEX]
                        align2_length = sum([int(i) for i in re.findall(r'\d+', CIGAR)])
                        align2_flag = int(content[FLAG_INDEX])
                        align2_pos = int(content[POS_INDEX])
                        # supplementary alignment
                        # align name; align chr; align length; align flag; same junction type; same geneid
                        if align1_name == align2_name and align1_chr == align2_chr and align1_length == align2_length and align2_flag == align1_flag+2048:
                            align2_aligned_seg = int(re.findall(r'\d+', CIGAR)[1])
                            align2_pos_plusM = align2_pos + align2_aligned_seg - 1
                            if align2_pos in junction_dict[align2_chr]['down']:
                                linear_align2 = True
                                align2_gene_id.extend(junction_dict[align2_chr]['down'][align2_pos]['gene_id'])
                            if align2_pos_plusM in junction_dict[align2_chr]['up']:
                                circ_align2 = True
                                align2_gene_id.extend(junction_dict[align2_chr]['up'][align2_pos_plusM]['gene_id'])
                            if (linear_align1 and linear_align2) or (circ_align1 and circ_align2):
                                align_gene_ids = list(set(align1_gene_id) & set(align2_gene_id))
                                if len(align_gene_ids) > 0: 
                                    for align_gene_id in align_gene_ids:
                                        ciri_result[align_gene_id] += 1
                else: # SM
                    if re.match(r'^\d+M\d+[SH]$', CIGAR) is not None:
                        align2_name = content[NAME_INDEX]
                        align2_chr = content[CHROM_INDEX]
                        align2_length = sum([int(i) for i in re.findall(r'\d+', CIGAR)])
                        align2_flag = int(content[FLAG_INDEX])
                        align2_pos = int(content[POS_INDEX])
                        if align1_name == align2_name and align1_chr == align2_chr and align1_length == align2_length and align2_flag == align1_flag+2048:
                            align2_aligned_seg = int(re.findall(r'\d+', CIGAR)[0])
                            align2_pos_plusM = align2_pos + align2_aligned_seg - 1
                            if align2_pos_plusM in junction_dict[align2_chr]['up']:
                                linear_align2 = True
                                align2_gene_id.extend(junction_dict[align2_chr]['up'][align2_pos_plusM]['gene_id'])
                            if align2_pos in junction_dict[align2_chr]['down']:
                                circ_align2 = True
                                align2_gene_id.extend(junction_dict[align2_chr]['down'][align2_pos]['gene_id'])
                            if (linear_align1 and linear_align1) or (circ_align1 and circ_align2):
                                align_gene_ids = list(set(align1_gene_id) & set(align2_gene_id))
                                if len(align_gene_ids) > 0: 
                                    for align_gene_id in align_gene_ids:
                                        ciri_result[align_gene_id] += 1
            else:
                align1_gene_id = []
                linear_align1 = False
                circ_align1 = False
                align_pattern = None
                align1 = ''
                if re.match(r'^\d+M\d+[SH]$', CIGAR):
                    align1_aligned_seg = int(re.findall(r'\d+', CIGAR)[0])
                    align1_pos = int(content[POS_INDEX])
                    align1_pos_plusM = align1_pos + align1_aligned_seg - 1
                    align1_chr = content[CHROM_INDEX]
                    if align1_pos_plusM not in junction_dict[align1_chr]['up'] and align1_pos in junction_dict[align1_chr]['down']:
                        continue
                    if align1_pos_plusM in junction_dict[align1_chr]['up']:
                        linear_align1 = True
                        for gene_id in junction_dict[align1_chr]['up'][align1_pos_plusM]['gene_id']:
                            if gene_id in ciri_result:
                                align1_gene_id.append(gene_id)
                    if align1_pos in junction_dict[align1_chr]['down']:
                        circ_align1 = True
                        for gene_id in junction_dict[align1_chr]['down'][align1_pos]['gene_id']:
                            if gene_id in ciri_result:
                                align1_gene_id.append(gene_id)
                    if len(align1_gene_id) == 0:
                        continue
                    first_align = True
                    align_pattern = 'MS'
                    align1_name = content[NAME_INDEX]
                    align1_length = sum([int(i) for i in re.findall(r'\d+', CIGAR)])
                    align1_flag = int(content[FLAG_INDEX])
                    align1 = line
                    continue
                if re.match(r'^\d+[SH]\d+M$', CIGAR):
                    align1_aligned_seg = int(re.findall(r'\d+', CIGAR)[1])
                    align1_pos = int(content[POS_INDEX])
                    align1_pos_plusM = align1_pos + align1_aligned_seg - 1
                    align1_chr = content[CHROM_INDEX]
                    if align1_pos not in junction_dict[align1_chr]['down'] and align1_pos_plusM in junction_dict[align1_chr]['up']:
                        continue
                    if align1_pos in junction_dict[align1_chr]['down']:
                        linear_align1 = True
                        for gene_id in junction_dict[align1_chr]['down'][align1_pos]['gene_id']:
                            if gene_id in ciri_result:
                                align1_gene_id.append(gene_id)
                    if align1_pos_plusM in junction_dict[align1_chr]['up']:
                        circ_align1 = True
                        for gene_id in junction_dict[align1_chr]['up'][align1_pos_plusM]['gene_id']:
                            if gene_id in ciri_result:
                                align1_gene_id.append(gene_id)
                    if len(align1_gene_id) == 0:
                        continue
                    first_align = True
                    align_pattern = 'SM'
                    align1_name = content[NAME_INDEX]
                    align1_length = sum([int(i) for i in re.findall(r'\d+', CIGAR)])
                    align1_flag = int(content[FLAG_INDEX])
                    align1 = line
                    continue
    
    return 

def main():
    use_list = sys.argv
    input_file = use_list[1]
    outputdir = use_list[2]
    gtffile = use_list[3]

    junction_dict = gtf_junction_extract(gtf_file)
    print('gtf reading finished.')
    input = pd.read_csv(input_file, sep='\t')
    
    for i in range(input.shape[0]):
        sample_id = input.iloc[i, 0]
        ciri_fn = input.iloc[i, 1]
        sam_fn = input.iloc[i, 3]
        ciri_gene_result, ciri_circ_result = read_ciri_result(ciri_fn)
        print('ciriquant result reading finished.')
        extract_juncreads_from_sam(sam_fn, ciri_gene_result, junction_dict, outputdir)
        output = '%s/%s.splicing_amount.output' % (outputdir, sample_id)
        with open(output, 'w') as f:
            f.write('circ_id\tgene_id\tsplicing_amount\n')
            for circ_id in ciri_circ_result:
                gene_ids = ciri_circ_result[circ_id]
                for gene_id in gene_ids:
                    f.write('{0}\t{1}\t{2}\n'.format(circ_id, gene_id, ciri_gene_result[gene_id]))
    return 

if __name__ == '__main__':
    main()
#!/usr/bin/env python

import argparse
import re

# python Ferina_deduper.py -f /input_sam_test.sam -o .sam  -u STL96.txt
# python Ferina_deduper.py -f /projects/bgmp/rferina/bioinfo/Bi624/dedup/Deduper-rferina/input_sam_test.sam -o deduplicated.sam  -u STL96.txt
def get_args():
    parser = argparse.ArgumentParser(description="Takes in a uniquely mapped SAM file, presorted by chromosome and position via samtools. Assumes reads are single-end.")
    parser.add_argument('-f', '--file', help='specify absolute file path of presorted input SAM file')
    parser.add_argument('-o', '--outfile', help='specify absolute file path of sorted output SAM file')
    parser.add_argument('-u', '--umi', help='specify file containing known valid unique molecular identifiers (UMIs)')
    return parser.parse_args()

args = get_args()

def strand_checker(bitwise):
    """
    Takes in SAM file bitwise flag. Returns strandedness (+ or -) if
    the read is mapped.
    """
    bitwise = int(bitwise)
    if ((bitwise & 16) == 16):
        rev_comp = True
        strandedness = 'negative'
    else:
        strandedness = 'positive' 
    return strandedness


# input: 107;0
# expected output: positive
# input: 83;16
# expected output: negative

def adjust_position(sam_position, strand, cigar):
    """
    Takes in the SAM file position, DNA strand, and CIGAR 
    string, parses the CIGAR string and returns the 5'
    start position.
    """
    # by default
    five_prime_start = int(sam_position)
    # if strand is positive, parse leftmost S if there
    if (strand == 'positive') & ('S' in cigar):
        left_s = re.findall(r"((\d+)[S])", cigar)[0]
        five_prime_start -= int(left_s[1])
    # if strand is negative, parse skipped region, deletion, and rightmost S
    elif strand == 'negative':
        if 'N' in cigar:
            skip = re.findall(r"((\d+)[N])", cigar)
            skipped_list = [int(i[1]) for i in skip]
            skipped = sum(skipped_list)
        else:
            skipped = 0
        if 'D' in cigar:
            deletion = re.findall(r"((\d+)[D])", cigar)
            deletion_list = [int(i[1]) for i in deletion]
            deleted = sum(deletion_list)
        else:
            deleted = 0
        if 'M' in cigar:
            match = re.findall(r"((\d+)[M])", cigar)
            match_list = [int(i[1]) for i in match]
            matched = sum(match_list)
        else:
            matched = 0
        if 'S' in cigar[-1]:
            right_s = re.findall(r"((\d+)[S]$)", cigar)[0]
            right_s = int(right_s[1])
        else:
            right_s = 0
        five_prime_start += skipped + deleted + right_s + matched
    return five_prime_start


# input: 5, +, 5S15M   
# expected output: 5

# input: 10, -, 3S17M
# expected output: 47

# input: 10, -, 4S10M50N5M3S
# expected output: 78

################################################################################################################################################################
# open files to write out to 
dedupe = open('deduplicated.sam', 'w')
dupes = open('duplicates.sam', 'w')
bad_umi = open('invalid_umis.sam', 'w')

# open umi file, place umis in set
with open(args.umi) as umi_file:
    valid_umis = set()
    for line in umi_file:
        line = line.strip('\n')
        valid_umis.add(line)


valid_reads = set()

with open(args.file) as input_sam:
    header_count = 0
    record_count = 0
    dup_count = 0
    unique_count = 0
    invalid_umi_count = 0
    for line in input_sam:
        # write out header lines
        if line.startswith('@'):
            dedupe.write(line)
            dupes.write(line)
            bad_umi.write(line)
            header_count += 1
        else:
            record_count += 1
            full_record = line.strip('\n')
            columns = full_record.split('\t')
            # extract relevant variables, the indexing starts at 0 which is inconsistent with columns of SAM file
            umi = columns[0].split(':')[-1]
            bitwise = columns[1]
            chromo = columns[2]
            sam_position = columns[3]
            cigar = columns[5]
            # check if umi is valid
            if umi not in valid_umis:
                # write out to incorrect umi sam file
                bad_umi.write(line)
                invalid_umi_count += 1
            # umi is valid
            else:
                # check if positive or negative strand
                strand = strand_checker(bitwise)
                # convert sam position to 5' start position
                five_prime_pos = adjust_position(sam_position, strand, cigar)
                # generate tuple with requirements 
                location_tuple = (umi, chromo, strand, five_prime_pos)
                # if tuple is not in the set, its unique; not a duplicate; add it to the set
                if location_tuple not in valid_reads:
                    # add to dict, increment counter
                    valid_reads.add(location_tuple)
                    # write out to deduplicated sam file
                    dedupe.write(line)
                    unique_count += 1
                # TEST: CHANGED IF TO ELIF
                # if tuple in set already, its a duplicate
                elif location_tuple in valid_reads:
                    # write out to duplicates file
                    dupes.write(line)
                    dup_count += 1
                   
    print('Header Lines:' , header_count)
    print('Total Records:', record_count)
    print('Unique Records:', unique_count)
    print('Duplicate Records:', dup_count)
    print('Invalid UMIs:', invalid_umi_count)

# close output files
dedupe.close()
dupes.close()
bad_umi.close()       
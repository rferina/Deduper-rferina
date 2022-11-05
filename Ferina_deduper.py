#!/usr/bin/env python

import argparse
import bioinfo
# FIX ARGPARSE FOR FILES

# python Ferina_deduper.py -f /projects/bgmp/rferina/bioinfo/Bi624/dedup/Deduper-rferina/test.sam -o large_test_output.sam  -u STL96.txt
# python Ferina_deduper.py -f /projects/bgmp/rferina/bioinfo/Bi624/dedup/Deduper-rferina/input_sam_test.sam -o deduplicated.sam  -u STL96.txt
def get_args():
    parser = argparse.ArgumentParser(description="Takes in a uniquely mapped SAM file, presorted by chromosome and position via samtools. Assumes reads are single-end.")
    parser.add_argument('-f', '--file', help='specify absolute file path of presorted input SAM file')
    parser.add_argument('-o', '--outfile', help='specify absolute file path of sorted output SAM file')
    parser.add_argument('-u', '--umi', help='specify file containing known valid unique molecular identifiers (UMIs)')
    return parser.parse_args()

args = get_args()

# open files to write out to 
dedupe = open('deduplicated.sam', 'w')
# dupes = open('duplicates.sam', 'w')
# bad_umi = open('invalid_umis.sam', 'w')

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
                # bad_umi.write(line)
                invalid_umi_count += 1
            # umi is valid
            else:
                # check if positive or negative strand
                strand = bioinfo.strand_checker(bitwise)
                # convert sam position to 5' start position
                five_prime_pos = bioinfo.adjust_position(sam_position, strand, cigar)
                # generate tuple with requirements 
                location_tuple = (umi, chromo, strand, five_prime_pos)
                # if tuple is not in the set, its unique; not a duplicate; add it to the set
                if location_tuple not in valid_reads:
                    # add to dict, increment counter
                    valid_reads.add(location_tuple)
                    # write out to deduplicated sam file
                    dedupe.write(line)
                    unique_count += 1
                # if tuple in set already, its a duplicate
                elif location_tuple in valid_reads:
                    # write out to duplicates file
                    # dupes.write(line)
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
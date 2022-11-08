#!/usr/bin/env python

import argparse
import bioinfo

def get_args():
    parser = argparse.ArgumentParser(description="Takes in a uniquely mapped SAM file, which can be presorted by chromosome and position via samtools. Assumes reads are single-end. Outputs deduplicated SAM file.")
    parser.add_argument('-f', '--file', help='specify absolute file path of presorted input SAM file')
    parser.add_argument('-o', '--outfile', help='specify absolute file path of sorted output SAM file')
    parser.add_argument('-u', '--umi', help='specify file containing known valid unique molecular identifiers (UMIs)')
    return parser.parse_args()

args = get_args()

# open file to write deduplicated SAM file out to 
dedupe = open(args.outfile, 'w')

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
            header_count += 1
        else:
            record_count += 1
            full_record = line.strip('\n')
            columns = full_record.split('\t')
            # extract relevant variables for determining duplicates, 
            # note the indexing starts at 0 which is inconsistent with columns of SAM file
            umi = columns[0].split(':')[-1]
            bitwise = columns[1]
            chromo = columns[2]
            sam_position = columns[3]
            cigar = columns[5]
            # check if umi is valid, increment invalid umi counter
            if umi not in valid_umis:
                invalid_umi_count += 1
            # now that umi is valid, continue checking if its a duplicate
            else:
                # check if positive or negative strand
                strand = bioinfo.strand_checker(bitwise)
                # convert sam position to 5' start position
                five_prime_pos = bioinfo.adjust_position(sam_position, strand, cigar)
                # generate tuple with requirements for duplicates
                location_tuple = (umi, chromo, strand, five_prime_pos)
                # if tuple is not in the set, its unique; not a duplicate
                if location_tuple not in valid_reads:
                    # add unique read to set, increment counter
                    valid_reads.add(location_tuple)
                    unique_count += 1
                    # write out full read to deduplicated sam file
                    dedupe.write(line)
                # if tuple in set already, its a duplicate, increment counter
                elif location_tuple in valid_reads:
                    dup_count += 1
                   
    print('Header Lines:' , header_count)
    print('Records Count:', record_count)
    print('Unique Records:', unique_count)
    print('Duplicate Records:', dup_count)
    print('Invalid UMIs:', invalid_umi_count)


dedupe.close()
       
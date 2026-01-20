#!/usr/bin/env python3

################################################################################
# MazF cleavage site enrichment analysis
# Copyright (C) 2026 New England Biolabs
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
################################################################################

import argparse
import pysam
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq')
    parser.add_argument('--output-file', type=str, default='read_count.csv')
    parser.add_argument('--output-fasta', type=str, default='trimmed_reads.fa')
    parser.add_argument('--label', type=str, default='count')
    parser.add_argument('--preset', type=str, default='CAA,7,TAA')
    args = parser.parse_args()

    # parse preset
    seq5, seqlen, seq3 = args.preset.split(',')
    seqlen = int(seqlen)

    counts = {}

    fasta_out = open(args.output_fasta, 'wt')

    with pysam.FastxFile(args.fastq) as fh:
        for entry in fh:
            if entry.sequence.find('N') != -1:
                continue

            p1 = -1
            p2 = -1
            offset = 0

            if seq5 != '':
                p1 = entry.sequence.find(seq5)

                if p1 == -1:
                    continue
                else:
                    p1 += len(seq5)
            else:
                p1 = 0

            if seq3 != '':
                p2 = entry.sequence.rfind(seq3)
                
                if p2 == -1:
                    continue
            else:
                p2 = len(entry.sequence)

            if p2 - p1 != seqlen:
                continue
            
            read = entry.sequence[p1:p2]

            if read not in counts:
                counts[read] = 0

            counts[read] += 1

            print('>' + entry.name + ' ' + entry.comment, file=fasta_out)
            print(read, file=fasta_out)

    fasta_out.close()

    # write read counts to output file
    f = open(args.output_file, 'wt')

    print('read', '%s' % args.label, sep=',', file=f)

    for read in sorted(counts.keys(), key=lambda x: counts[x]):
        print(read, counts[read], sep=',', file=f)

    f.close()

if __name__ == '__main__':
    sys.exit(main())

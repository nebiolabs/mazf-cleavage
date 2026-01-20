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
import pandas as pd
import numpy as np
import random
import sys
import timeit

def calculate_enrichment(control_lib, sample_lib):
    ### merge by read sequence
    f = control_lib.merge(sample_lib, on='read', how='outer')
    f.columns = ['read','control','sample']

    ### if some reads are not present in control/sample init with 0
    f.fillna(0, inplace=True)

    ### add pseudocount (+1)
    f = f.astype({'control':int, 'sample':int})
    f['control'] = f['control'] + 1
    f['sample']  = f['sample']  + 1

    ### normalize read counts in control and sample
    f['control_norm'] = f['control'] / f['control'].sum()
    f['sample_norm']  = f['sample']  / f['sample'].sum()

    ### compute enrichment ratio
    f['ratio'] = f['control_norm'] / f['sample_norm']
    f['logratio'] = np.log10(f['ratio'])

    ### sort reads by enrichment ratio
    f.sort_values(by=['ratio'], ascending=False, inplace=True)

    return(f)

def sample(control_lib, sample_lib, max_it):
    # determine number of reads in the control and sample libraries
    control_lib_size = control_lib['count'].sum()
    sample_lib_size = sample_lib['count'].sum()

    # determine how many reads to sample
    n = sample_lib_size if sample_lib_size <= control_lib_size else control_lib_size

    print()
    print('Control library           : %9i' % control_lib_size)
    print('Sample library            : %9i' % sample_lib_size)
    print('Bootstrap control library : %9i' % n)

    print()
    print('Bootstrap:')

    logratio = []

    for it in range(max_it):
        start = timeit.default_timer()

        # bootstrap
        rnd = random.choices(control_lib['read'], weights = control_lib['count'], k = n)

        # reformat bootstrap control library
        counts = {}

        for i,r in enumerate(rnd):
            if r not in counts:
                counts[r] = 0
            counts[r] += 1
            
        data = []

        for r in sorted(counts):
            data.append([r,counts[r]])

        control_lib_bootstrap = pd.DataFrame(data = data, columns = ['read','count'])

        f = calculate_enrichment(control_lib, control_lib_bootstrap)

        logratio.append(f.iloc[0]['logratio'])

        end = timeit.default_timer()

        print('iteration=%i max_logratio=%.7f time=%.f' % (it, f.iloc[0]['logratio'], end-start))
    
    return(logratio)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('control')
    parser.add_argument('sample')
    parser.add_argument('--bootstrap', type=int, default=100)
    parser.add_argument('--quantile', type=float, default=0.95)
    parser.add_argument('--output-csv', default='output.csv')
    parser.add_argument('--gen-fasta', action='store_true', default=False)
    parser.add_argument('--output-fasta', default='output.fasta')
    args = parser.parse_args()

    fc = pd.read_csv(args.control)
    fs = pd.read_csv(args.sample)

    ### load read count for control and actual sample
    fc.columns = ['read','count']
    fs.columns = ['read','count']

    logratio = sample(fc, fs, args.bootstrap)

    print()
    print('q95_logratio=%.7f' % np.quantile(a=logratio, q=0.95))
    print('q99_logratio=%.7f' % np.quantile(a=logratio, q=0.99))
    print('max_logratio=%.7f' % np.max(logratio))

    threshold = np.quantile(a=logratio, q=args.quantile)

    print()
    print('treshold=%.7f' % threshold)

    f = calculate_enrichment(fc, fs)

    ### save control/sample info to output file
    f.to_csv(args.output_csv, sep=',', index=False)

    if args.gen_fasta:
        counter = 1

        fout = open(args.output_fasta, 'wt')

        for index, row in f.iterrows():
            if row['logratio'] > threshold:
                print('enriched_sequence=%s ratio=%.7f logratio=%.7f' % (row['read'], row['ratio'], row['logratio']))
                
                for c in range(int(row['ratio'] * 10)):
                    print('>r' + str(counter), file=fout)
                    print(row['read']+'N', file=fout)
                    counter += 1

                    if counter > 500000:
                        break

            if counter > 500000:
                break

        fout.close()

if __name__ == '__main__':
    sys.exit(main())

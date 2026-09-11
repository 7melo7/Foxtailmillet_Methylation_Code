#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse
import glob
import gzip

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf
    '''
    print message


def mergefiles(dfs=[], on=''):
    """Merge a list of files based on one column"""
    if len(dfs) == 1:
         return "List only have one element."

    elif len(dfs) == 2:
        df1 = dfs[0]
        df2 = dfs[1]
        df = df1.merge(df2, on=on)
        return df

    # Merge the first and second datafranes into new dataframe
    df1 = dfs[0]
    df2 = dfs[1]
    df = dfs[0].merge(dfs[1], on=on)

    # Create new list with merged dataframe
    dfl = []
    dfl.append(df)

    # Join lists
    dfl = dfl + dfs[2:] 
    dfm = mergefiles(dfl, on)
    return dfm

def read_fasta(infile):
    fasta_seq = defaultdict(lambda : str())
    for record in SeqIO.parse(infile, "fasta"):
        fasta_seq[record.id] = str(record.seq)
        #print record.id, fasta_seq[record.id]
    return fasta_seq


#scaffold_1      4       +       0       0       CHH     CAC     1       1       NA
#scaffold_1      6       +       0       0       CHG     CAG     1       1       NA
def read_CX_report(infile, prefix):
    ofile_CG   = open('%s.CG.bedgraph' %(prefix), 'w')
    ofile_CHG  = open('%s.CHG.bedgraph' %(prefix), 'w')
    ofile_CHH  = open('%s.CHH.bedgraph' %(prefix), 'w')
    #print >> ofile_CG, 'track type=bedGraph'
    #print >> ofile_CHG, 'track type=bedGraph'
    #print >> ofile_CHH, 'track type=bedGraph'
    with gzip.open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if unit[9] == 'NA':
                    continue
                rate = float(unit[3])/(float(unit[3])+float(unit[4]))
                start = unit[1]
                end   = str(int(unit[1]) + 1)
                if unit[5] == 'CG':
                    print >> ofile_CG, '%s\t%s\t%s\t%s' %(unit[0], start, end, str(rate))
                elif unit[5] == 'CHG':
                    print >> ofile_CHG, '%s\t%s\t%s\t%s' %(unit[0], start, end, str(rate))
                elif unit[5] == 'CHH':
                    print >> ofile_CHH, '%s\t%s\t%s\t%s' %(unit[0], start, end, str(rate))
    ofile_CG.close()
    ofile_CHG.close()
    ofile_CHH.close()

#Cell    Subclone
#FEL024_S_AACGGGATCATGTCAG       Cluster1
#FEL024_S_AAAGTGACATGCCGAC       Cluster1
#FEL024_S_ACTTCGCGTATACCCA       Cluster4
def read_cell_anno(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    file_df.columns = ["Cell", "Anno"]
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        data[file_df['Cell'][i]] = file_df['Anno'][i]
    return data

def read_large_matrix(infile):
    # determine and optimize dtype
    # Sample 100 rows of data to determine dtypes.
    file_test = pd.read_csv(infile, sep="\t", header=0, nrows=100)
    float_cols = [c for c in file_test if file_test[c].dtype == "float64"]
    int_cols = [c for c in file_test if file_test[c].dtype == "int64"]
    if float_cols > 0:
        dtype_cols = {c: np.float32 for c in float_cols}
    elif int_cols > 0:
        dtype_cols = {c: np.int32 for c in int_cols}
    file_df = pd.read_csv(infile, sep="\t", header=0, dtype=dtype_cols)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    read_CX_report(args.input, args.output)

if __name__ == '__main__':
    main()


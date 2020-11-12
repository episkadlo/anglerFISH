#!/usr/bin/env python
scriptName = 'customBed2Fasta'
Version = 1.2


import pandas
import os
import argparse

def main():

    #loop over the rows of input bed file in order to extract the sequence alone from the seq column and save it to the new file file
    def bedToFasta(file, out):
            df = pandas.read_csv(str(file), sep = '\t', header = None)
            df.columns= ['chr', 'pos1', 'pos2', 'seq', 'Tm']
            probe_inx = []
            for i in range (1, len(df)+1): probe_inx.append('>probe_' + str(i) + "_" + str(out))
            new_inx = pandas.Series(probe_inx)
            out_df = pandas.concat([new_inx, df['seq']], axis = 1, sort = False)
            f= open(str(out) + '.fasta','w+')
            for index, row in out_df.iterrows():
                f.write(str(row[0])+"\n")
                f.write(str(row['seq'])+"\n")
            f.close()


    #Initialize parse_args to get input file
    parser = argparse.ArgumentParser(description="Extracts sequences from .bed file and saves them as a single .fasta file")

    # add ArgumentParser
    parser.add_argument("-f", "--File", required=True, help = "input bed file, from which the single sequences will be extracted and saved as a fasta file; reqired")
    parser.add_argument("-o", "--out", required=True, help = "output fasta file name; required")
    args = parser.parse_args()
    inFile = args.File
    outSuffix = args.out


    #call the bedToFasta on the user's input file to extract sequences
    bedToFasta(inFile, outSuffix)

if __name__ == '__main__':
    main()

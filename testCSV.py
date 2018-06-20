#!/usr/bin/env python
#-*- coding : utfs -*-

"""
Author: Pierre Merckaert
Contact: merckaert.pierre@gmail.com
Date: 06/2017

Description: Calculates the sgRNA efficiency score for a given 30-mer sgRNA or a csv file containing all sgRNA to predict.
Input: 30mer sgRNA+context sequence, NNNN[sgRNA sequence]NGGNNN
Output: efficiency score
"""


import csv, argparse, sys, subprocess, pickle, StandardScaler, VarianceThreshold
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

def csvParser(processed_sgRNAs):
    with open(processed_sgRNAs) as f:
        file = csv.reader(f,delimiter=',')
        gRNAs = []
        for row in file:
            gRNA = row
            gRNAs.append(gRNA)
    return(gRNAs)

def get_parser():
    parser = argparse.ArgumentParser(description='''Calculates the sgRNA efficiency score for a given 30-mer sgRNA OR a csv file containing all sgRNA to predict.''',
    epilog="")
    parser.add_argument('--seq',
        type=str,
        help='30mer sgRNA + context sequence, NNNN[sgRNA sequence]NGGNNN')
    parser.add_argument('--csv', type=argparse.FileType('r'),
        help='''csv file containing all sgRNA to predict under a column named "gRNA_30mer".\n
        Format : Comma-delimited csv file with the list of 30mer sgRNAs in the first column. \n
        A header row is required.''')
    return(parser)

if __name__ == '__main__':
    args = get_parser().parse_args()
    assert ((args.csv is not None and args.seq is None) or (args.csv is None and args.seq is not None)), "you have to specify either 30mer sgRNA sequence or a csv file (see help)"
    if args.csv != None : 

        # Define command and arguments
        command = 'Rscript'
        path2script = 'Rpreprocessing/sgRNA_Input_Processing.R'

        # Build subprocess command
        cmd = [command, path2script, args.csv.name]
        print(cmd)

        # check_output will run the command and store to result
        subprocess.run(cmd)

        df_gRNAs = np.array(csvParser('Processed_Sequence.csv'), dtype='<U32')
        print(df_gRNAs)
        names = df_gRNAs[0,:]
        gRNAs_param = np.array(df_gRNAs[1:,1:], dtype='float64')
        print(gRNAs_param)

        selector = VarianceThreshold()
        selector.fit_transform(gRNAs_param)
        names = names[selector.get_support(False)]
        gRNAs_param = gRNAs_param[:,selector.get_support(False)]
        print(gRNAs_param.shape)
        print(names.shape)

        sc = StandardScaler()
        gRNAs_param[:,0:24] = sc.fit_transform(gRNAs_param[:,0:24])


        model_file = '201FS_Ridge_model.pickle'

        try:
            with open(model_file, 'rb') as f:
                model,indexes = pickle.load(f)
        except:
            raise Exception("could not find model stored to file %s" % model_file)

        gRNAs_param = gRNAs_param[:,indexes]

        scores = model.predict(gRNAs_param)[0]
        print(scores)

    elif len(args.seq)==30 :
        gRNA = args.seq
        print(gRNA)
    else :
        print("sgRNA must be 30mer long (see help)")


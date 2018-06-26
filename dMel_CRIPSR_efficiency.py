#!/usr/bin/env python
#-*- coding : utfs -*-

"""
Author: Pierre Merckaert
Contact: merckaert.pierre@gmail.com
Date: 06/2017

Description: Calculates the sgRNA efficiency score for a given 30-mer sgRNA or a csv file containing all sgRNA to predict.
Input: 30mer sequence(s) : 20mer sgRNA + context sequence, NNNN[20mer sgRNA sequence]NGGNNN
Output: efficiency score(s)
"""

import csv, argparse, os, sys, subprocess, pickle, sklearn, re
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np

# Description: get the input arguments and sets up the --help
# Output: parsed arguments
def get_parser():
	parser = argparse.ArgumentParser(description='''Calculates the sgRNA efficiency score for a given 30-mer sgRNA OR a csv file containing all sgRNA to predict. The score ranges between 0 and 1. The higher the better the sgRNA efficiency is predicted to be.''', epilog="")
	parser.add_argument('--seq',
		type=str,
		help='30mer sgRNA + context sequence, NNNN[sgRNA sequence]NGGNNN')
	parser.add_argument('--csv', type=argparse.FileType('r'),
		help='''csv file containing all sgRNA to predict under a column header.\n
		Format : Comma-delimited csv file with the list of 30mer sgRNAs in the first column. \n
		A header row is required.''')
	parser.add_argument('--out', type=str,
		help='''Path and/or Name of the predictions output csv file. Default is PATH_TO_FOLDER/sgRNA_efficiency_prediciton/sgRNA_predictions.csv''')
	return(parser)

# Description: Run the R pipeline to extract the features from the input sequences
# Input: f_input is the csv file or the 30mer sequence 		path is the path of the current script
# Output: outputs the extracted features in /Rpreprocessing/R_Featurized_sgRNA.csv
def Rpreprocessing(f_input, path):
	command = 'Rscript'
	path2script = str(path)+'/Rpreprocessing/sgRNA_Input_Processing.R'
	cmd = [command, path2script, f_input]
	subprocess.run(cmd)

# Description: open a given csv file and appends its rows to an array
# Input: csv file of the extracted features from the R pipeline
# Output: array of the extracted features from the R pipeline
def csvParser(processed_sgRNAs):
	try:
		with open(processed_sgRNAs) as f:
			file = csv.reader(f,delimiter=',')
			gRNAs = []
			for row in file:
				gRNA = row
				gRNAs.append(gRNA)
	except:
		raise Exception("could not find input stored to file %s" % processed_sgRNAs)

	return(gRNAs)

# Description: Python pipeline to load the ML model and create the dataset for the efficiency prediction  
# Input: path of the current file
# Output: array of the 23 and 30mer sgRNAs, features of these sgRNAs and the ML model to make the efficiency predictions
def PythonPreprocessing(path):
	# Opens the extracted features of the R Pipeline
	# The file format is : 30mer(1), 23mer(1), GC%(1), Tm(4), indep Order1(4), indep Order2(16), dep Order1, dep Order2, NGGN
	df_gRNAs = np.array(csvParser(str(path)+'/Rpreprocessing/R_Featurized_sgRNA.csv'), dtype='<U32')

	#First line is the features' names
	names = df_gRNAs[0,2:]
	#Features of the sequences
	gRNAs_param = np.array(df_gRNAs[1:,2:], dtype='float64')

	#Standard Scaling of the non-binary features
	sc = StandardScaler()
	gRNAs_param[:,0:24] = sc.fit_transform(gRNAs_param[:,0:24])

	#Load the Ridge regression model and the index of relevant features (201 out of 625)
	model_file = str(path)+'/201FS_Ridge_model.pickle'
	try:
		pickle_in = open(model_file,"rb")
		p_load = pickle.load(pickle_in)
		RidgeModel = p_load['model']
		idx = p_load['df_indexes']
	except:
	    raise Exception("could not find model stored to file %s" % model_file)
	
	gRNAs_param = gRNAs_param[:,idx] #Keep only the 201 relevant features
	gRNAs_seq = df_gRNAs[1:,0:2] 	 #23 and 30mer sequences

	return(gRNAs_seq, gRNAs_param, RidgeModel)


if __name__ == '__main__':
	args = get_parser().parse_args()
	path = os.path.dirname(os.path.realpath(__file__)) #path of the directory of this python script

	assert ((args.csv is not None and args.seq is None) or (args.csv is None and args.seq is not None)), "you have to specify either 30mer sgRNA sequence OR a csv file (see help)"
	if args.csv != None : 
		# Extract features from input csv file
		Rpreprocessing(os.path.realpath(args.csv.name), path)
	elif (len(args.seq)==30) & (re.match(r"^[ATGCatgc]*$", args.seq)!=None):
		# Extract features from input sequence
		Rpreprocessing(args.seq.upper(), path)
		print_res = True #print prediction score in the terminal
	else :
		print("wrong input (see --help)")

	#get the sgRNAs sequences, their features and the Machine Learning model
	gRNAs_seq, gRNAs_param, RidgeModel = PythonPreprocessing(path)

	scores = RidgeModel.predict(gRNAs_param) #Predictions !

	#Change the output if specified
	if args.out != None :
		output = args.out
	else : 
		output = path+"/sgRNA_predictions.csv"

	#outputs the scores in a dataframe and save 
	pd.DataFrame({'gRNA_23mer' : gRNAs_seq[:,1], 'gRNA_30mer' : gRNAs_seq[:,0], 'scores' : scores}).to_csv(output, index=False)

	print('DONE. Results exported to %s' % output)

	#print prediction score in the terminal
	if print_res :
		print('Predicted efficiency score [0:1] for', args.seq.upper(),'= ',scores, '(the higher the better)')








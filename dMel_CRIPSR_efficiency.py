#!/usr/bin/env python
#-*- coding : utfs -*-

"""
Author: Pierre Merckaert
Contact: merckaert.pierre@gmail.com
Date: 06/2018

Description: Calculates the sgRNA efficiency score for a given 23 or 30-mer sgRNA or a csv file containing all sgRNA to predict.
Input: 23 or 30mer sequence(s) : 20mer sgRNA + PAM (+ context sequence), NNNN[20mer sgRNA sequence]NGGNNN
Output: efficiency score(s)
"""

import os
import utils.util as util

if __name__ == '__main__':
	args = util.get_parser().parse_args()
	path = os.path.dirname(os.path.realpath(__file__)) #path of the directory of this python script

	#Launch the R preprocessing pipeline depending of the input arguments
	print_res = util.Launch_R_Preprocessing(args,path)

	#launch the Python pipeline to retrieve the sgRNAs sequences, their features and the Machine Learning model
	gRNAs_seq, gRNAs_param, ML_Model = util.PythonPreprocessing(path)

	#Predictions!
	scores = ML_Model.predict(gRNAs_param) 

	#Output and save results depending on arguments input and processing outputs
	util.Output_Results(gRNAs_param, gRNAs_seq, scores, args, path)

	#print prediction score in the terminal
	if print_res :
		print('Predicted efficiency score [0:1] for', args.seq.upper(),'= ',scores, '(the higher the better)')

#!/usr/bin/env python
#-*- coding : utfs -*-

"""
Author: Pierre Merckaert
Contact: merckaert.pierre@gmail.com
Date: 06/2018

Description: Calculates the sgRNA efficiency score for a given 20, 23 or 30-mer sgRNA or a csv file containing all sgRNA to predict.
Input: 20, 23 or 30mer sequence(s) : 20mer sgRNA + PAM (+ context sequence), NNNN[20mer sgRNA sequence]NGGNNN
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
	gRNAs_seq, gRNAs_param_REG, gRNAs_param_CLASS, ML_Model_REG, ML_Model_CLASS = util.PythonPreprocessing(path)
	
	#Predictions!
	scores_REG = ML_Model_REG.predict(gRNAs_param_REG)
	scores_CLASS = ML_Model_CLASS.predict(gRNAs_param_CLASS)

	#Output and save results depending on arguments input and processing outputs
	util.Output_Results(gRNAs_seq, gRNAs_param_REG, scores_REG, gRNAs_param_CLASS, scores_CLASS, args, path)

	#print prediction score in the terminal
	if print_res :
		if args.bin != "no":
			print('Predicted efficiency score [0:1] for', args.seq.upper(),'= ',scores_REG, '(the higher the better)')
			if scores_CLASS == 0:
				print('The classification model predicts that ', args.seq.upper(),' will NOT be an effective sgRNA.')
			else : print('The classification model predicts that ', args.seq.upper(),' will be an effective sgRNA.')
		else : print('Predicted efficiency score [0:1] for', args.seq.upper(),'= ', scores_REG, '(the higher the better)')

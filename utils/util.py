#!/usr/bin/env python
#-*- coding : utfs -*-

"""
Author: Pierre Merckaert
Contact: merckaert.pierre@gmail.com
Date: 06/2018
"""

import csv, argparse, os, sys, subprocess, pickle, sklearn, re
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from collections import OrderedDict

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
		help='''Path and/or Name of the predictions output csv file. \n
		Default is PATH_TO_FOLDER/sgRNA_efficiency_prediciton/sgRNA_predictions.csv''')
	parser.add_argument('--bin', type=str, default='no',
		help='''Choose with "yes" or "no" (default "no") whether or not to include binary classification prediction of the sgRNA efficiency in the output csv file''')
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
	try :
		df_gRNAs = np.array(csvParser(str(path)+'/Rpreprocessing/R_Featurized_sgRNA.csv'), dtype='<U32')
	except :
		raise Exception("could not find Featurized sgRNAs stored to file %s" % df_gRNAs)	

	if df_gRNAs.shape[1] == 470 : #number of features extracted from 23mer sgRNAs
		#Features of the sequences
		gRNAs_param = np.array(df_gRNAs[1:,1:], dtype='float64')
		#Load the Ridge regression model and the index of relevant features
		model_REG_file = str(path)+'/utils/23mer_135FS_Ridge_REGmodel.pickle'
		model_CLASS_file = str(path)+'/utils/23mer_120FS_GBC_BINmodel.pickle'

		gRNAs_seq = df_gRNAs[1:,0] #23mer sequences

	elif df_gRNAs.shape[1] == 410 : #number of features extracted from 20mer sgRNAs
		#Features of the sequences
		gRNAs_param = np.array(df_gRNAs[1:,1:], dtype='float64')
		#Load the Ridge regression model and the index of relevant features
		model_REG_file = str(path)+'/utils/20mer_326FS_GBRT_REGmodel.pickle'
		model_CLASS_file = str(path)+'/utils/20mer_137FS_GBC_BINmodel.pickle'

		gRNAs_seq = df_gRNAs[1:,0] #20mer sequences

	elif df_gRNAs.shape[1] == 627 :	#number of features extracted from 30mer sgRNAs
		#Features of the sequences
		gRNAs_param = np.array(df_gRNAs[1:,2:], dtype='float64')
		#Load the Ridge regression model and the index of relevant features
		model_REG_file = str(path)+'/utils/30mer_121FS_LinSVR_REGmodel.pickle'
		model_CLASS_file = str(path)+'/utils/30mer_400FS_GBRC_BINmodel.pickle'

		gRNAs_seq = df_gRNAs[1:,0:2] 	 #23 and 30mer sequences
	else :
		raise Exception('Error : wrong R pipeline output dimension')

	#Standard Scaling of the non-binary features
	sc = StandardScaler()
	gRNAs_param[:,0:25] = sc.fit_transform(gRNAs_param[:,0:25])

	#Load regression model
	try:
		pickle_in_REG = open(model_REG_file,"rb")
		p_load_REG = pickle.load(pickle_in_REG)
		ML_Model_REG = p_load_REG['model']
		idx_REG = p_load_REG['df_indexes']
	except:
	    raise Exception("could not find Regression model stored to file %s" % model_REG_file)

	#Load classification model
	try:
		pickle_in_CLASS = open(model_CLASS_file,"rb")
		p_load_CLASS = pickle.load(pickle_in_CLASS)
		ML_Model_CLASS = p_load_CLASS['model']
		idx_CLASS = p_load_CLASS['df_indexes']
	except:
	    raise Exception("could not find classification model stored to file %s" % model_CLASS_file)

	gRNAs_param_REG = gRNAs_param[:,idx_REG] #Keep only the relevant features
	gRNAs_param_CLASS = gRNAs_param[:,idx_CLASS]

	return(gRNAs_seq, gRNAs_param_REG, gRNAs_param_CLASS, ML_Model_REG, ML_Model_CLASS)


# Description: run the R pipeline depending on the input given 
# Input: arguments given to the script and path of the current file
# Output: csv of the featurized sgRNAs and if the prediction scores should be written in the terminal
def Launch_R_Preprocessing(args,path):
	assert ((args.csv is not None and args.seq is None) or (args.csv is None and args.seq is not None)), "you have to specify either 30mer sgRNA sequence OR a csv file (see help)"
	if args.csv != None : 
		# Extract features from input csv file
		Rpreprocessing(os.path.realpath(args.csv.name), path)
		print_res = False
	elif ((len(args.seq)==30) | (len(args.seq)==23) | (len(args.seq)==20)) & (re.match(r"^[ATGCatgc]*$", args.seq)!=None):
		# Extract features from input sequence
		Rpreprocessing(args.seq.upper(), path)
		print_res = True #print prediction score in the terminal
	else :
		raise Exception("wrong input (see --help)")

	return(print_res)


# Description: Outputs the sgRNAs prediction scores in csv files
# Input: 
	#gRNAs_param_REG : values of the relevant features for regression model
	#gRNAs_seq : list of the sgRNAs to predict
	#scores : prediction scores
# Output: save the prediction scores to csv file
def Output_Results(gRNAs_seq, gRNAs_param_REG, scores_REG, gRNAs_param_CLASS, scores_CLASS, args, path):
	if gRNAs_param_REG.shape[1] == 135 : #shape of the dataframe after feature selection for 23mer sgRNAs
		#Change the output if specified
		if args.out != None :
			output = args.out
		else : 
			output = path+"/23mer_sgRNA_predictions.csv"
		#outputs the scores in a dataframe and save
		try :
			if args.bin == "no":
				pd.DataFrame(OrderedDict({'gRNA_23mer' : gRNAs_seq, 'regression scores' : scores_REG})).to_csv(output, index=False)
			else:
				pd.DataFrame(OrderedDict({'gRNA_23mer' : gRNAs_seq, 'regression scores' : scores_REG, 'classification scores' : scores_CLASS})).to_csv(output, index=False)
		except:
			raise Exception("Error in saving prediction results to %s. CSV file needed." % output)

	elif gRNAs_param_REG.shape[1] == 326 : #shape of the dataframe after feature selection for 23mer sgRNAs
		#Change the output if specified
		if args.out != None :
			output = args.out
		else : 
			output = path+"/20mer_sgRNA_predictions.csv"
		#outputs the scores in a dataframe and save 
		try :
			if args.bin == "no":
				pd.DataFrame(OrderedDict({'gRNA_20mer' : gRNAs_seq, 'regression scores' : scores_REG})).to_csv(output, index=False)
			else:	
				pd.DataFrame(OrderedDict({'gRNA_20mer' : gRNAs_seq, 'regression scores' : scores_REG, 'classification scores' : scores_CLASS})).to_csv(output, index=False)
		except:
			raise Exception("Error in saving prediction results to %s. CSV file needed." % output)
		
	elif gRNAs_param_REG.shape[1] == 121 : #shape of the dataframe after feature selection for 30mer sgRNAs
		#Change the output if specified
		if args.out != None :
			output = args.out
		else : 
			output = path+"/30mer_sgRNA_predictions.csv"
		#outputs the scores in a dataframe and save 
		try :
			if args.bin == "no":
				pd.DataFrame(OrderedDict({'gRNA_23mer' : gRNAs_seq[:,1], 'gRNA_30mer' : gRNAs_seq[:,0], 'regression scores' : scores_REG})).to_csv(output, index=False)
			else:	
				pd.DataFrame(OrderedDict({'gRNA_23mer' : gRNAs_seq[:,1], 'gRNA_30mer' : gRNAs_seq[:,0], 'regression scores' : scores_REG, 'classification scores' : scores_CLASS})).to_csv(output, index=False)	
		except:
			raise Exception("Error in saving prediction results to %s. CSV file needed." % output)
	else :
		raise Exception('Error : wrong R pipeline output dimension')
	print('DONE. Results exported to %s' % output)
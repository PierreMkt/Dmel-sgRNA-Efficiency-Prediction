# sgRNA-Efficiency-Prediction

Python and R pipeline to predict sgRNA efficiency via a machine learning algorithm (Ridge regression) trained and optimized on CRIPSR pooled gene essentiality screens.   
The score returned ranges between 0 and 1, the higher the more effective the sgRNA is predicted to be.

### Requirements
Softwares : Python3 & R   
Python Packages :
   1. sklearn
   2. pickle    
   3. pandas    
   4. numpy
   5. scipy

## Inputs
1/ _--seq_ : the raw 30mer sequence.   
Format of the raw sequence :  
20mer sgRNA + PAM + context sequence : **NNNN[20mer sgRNA sequence]NGGNNN**   
      
2/ _--csv_ : .csv file with a header and the 30mer sgRNAs beneath   
Format of the **Comma delimited** .csv file

|  _Header_  |
| ------------- |
|  GGTTGCAGCTTTAGTGGTCGACAACGGATC  |
|  AAATCGCGACCGGCTAGATCCAAACGGAGG  | 
|              ...                 | 
|  GCCCAACCATGGGCAAGCGTCTGCAGGACA  | 

3/ _--out_ : (optional argument) Path and/or Name of the predictions output csv file. Default is PATH_TO_FOLDER/sgRNA_efficiency_prediciton/sgRNA_predictions.csv

## Outputs
The R Pipeline outputs the dataframe of the extracted features from the 30mer sequences in PATH_TO_FOLDER/sgRNA_efficiency_prediciton/Rpreprocessing/R_Featurized_sgRNA.csv

The prediction score for each sgRNA is written in a csv file that can be specified via the _--out_ argument. By default, the scores are exported to PATH_TO_FOLDER/sgRNA_efficiency_prediciton/sgRNA_predictions.csv

## Running the code
Preferably, make the "sgRNA_efficiency_prediciton" folder your working directory.

To run model prediction on a raw sgRNA, type on the terminal:   
_python dMel_CRIPSR_efficiency.py --seq <30mer sequence> [--out <PATH_TO_OUTPUT>]_
   
To run model prediction on a csv file, type on the terminal:    
_python dMel_CRIPSR_efficiency.py --csv <PATH_TO_CSV_FILE> [--out <PATH_TO_OUTPUT>]_

### Examples
_python dMel_CRIPSR_efficiency.py --seq TGGAGGCTGCTTTACCCGCTGTGGGGGCGC_    
Output : _Predicted efficiency score [0:1] for TGGAGGCTGCTTTACCCGCTGTGGGGGCGC =  [0.38789814] (the higher the better)_

_python dMel_CRIPSR_efficiency.py --csv sgRNA_to_predict.csv_     
Output : _DONE. Results exported to /home/.../sgRNA_efficiency_prediction/sgRNA_predictions.csv_

## Authors

**Pierre Merckaert** - [PierreMkt](https://github.com/PierreMkt)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments



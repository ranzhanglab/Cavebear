# Cross-species pseudotime prediction

<img width="2210" height="772" alt="Cavebear_overview_github" src="https://github.com/user-attachments/assets/05f8c333-06b7-4750-98eb-853acab1667d" />  

Concept: Cavebear utilizes an accurately labled, densley colleceted scRNA-seq reference dataset to predict the pseudotime of a target scRNA-seq dataset that has incomplete or inaccurate time-labels. The reference and target datasets can be from different species and/or experimental systems (i.e. _in vivo_ and _in vitro_).

Method: Cavebear aligns species or experimental systems (step 1) and uses the learned cell embeddings of the reference cells to train the time predictor (step 2). The trained time predictor is subsequently applied to the learned cell embeddings of the target species or system from step 1 to predict the pseudotime of the target cells (step 3). 

## Installation:
Install through conda:
```
conda env create -f environment.yml
conda activate
```

## Example Run:
```
cd ./src/
bash ./run.sh
```

## Basic Usage:
### 1. Cross-species alignment
```
python ./src/cavebear_pytorch_cvae.py --input_h5ad {path/to/input_file.h5ad} --predict train
```
example input file: ./data/example.h5ad
For hyperparameter tuning, recommend tuning the learning_rate (--learining-rate {0.01 0.001 0.0001}). Other potential hyperparameters to test are nlayers, embed_dim.  
> Defaults:  
&nbsp;&nbsp;--learning_rate 0.001  
&nbsp;&nbsp;--nlayers 3  
&nbsp;&nbsp;--embed_dim 25  

When predicting pseudotime across experimental systems (i.e. __in vitro_ and _in vivo_), include the discriminator argument and tune the discriminator weight.
> Arguments to include for discriminator:  
&nbsp;&nbsp;--dis dis   
&nbsp;&nbsp;--discriminator_weight {1.0 2.0 5.0 10.0} 


### 2. Model selection

### 3. Psuedotime training and prediction
```
python /src/cavebear_pytorch_cvae.py --input_h5ad /data/example.h5ad --predict predict --train_species mouse --target_species zebrafish --time_label mouse_age {optional: arguments for hyperparameters of best model if not default}
```
Where  `time_label` is the name of the column containing the training species time data (ie mouse age when cells were collected).  
If the best model from cross-species alignment has hyperparameters that differ from the defaults, you must include them as arguments.  


### 3. (OPTIONAL) Extract gene probabilites after cross-species alignment
```
python /src/cavebear_pytorch_cvae.py --input_h5ad /data/example.h5ad --predict px_decoder {optional: arguments for hyperparameters of best model if not default}
```

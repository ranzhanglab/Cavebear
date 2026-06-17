# Cross-species pseudotime prediction

<img width="2210" height="772" alt="Cavebear_overview_github" src="https://github.com/user-attachments/assets/05f8c333-06b7-4750-98eb-853acab1667d" />  

Concept: Cavebear utilizes an accurately labled, densley colleceted scRNA-seq reference dataset to predict the pseudotime of a target scRNA-seq dataset that has incomplete or inaccurate time-labels. The reference and target datasets can be from different species and/or experimental systems (i.e. _in vivo_ and _in vitro_).

Method: Cavebear aligns species or experimental systems (step 1) and uses the learned cell embeddings of the reference cells to train the time predictor (step 2). The trained time predictor is subsequently applied to the learned cell embeddings of the target species or system from step 1 to predict the pseudotime of the target cells (step 3). 

## Installation:
Install through conda:
```
conda env create -f environment.yml
conda activate cavebear
```

## Example Run:
```
bash ./src/run.sh --input {input.h5ad} --train_species {species1} --target_species {species2} [OPTIONAL ARGUMENTS]
```

## Input/ Preprocessing
Required input is a scRNA-seq (*.h5ad format) containing a combined sample X gene matrix from a reference dataset and a target dataset. The h5ad file will be read in as an anndata and the metadata will be located in adata.obs as a pandas dataframe. Metadata to include are at minimum 'species', 'batch', and 'age' or 'time', with optional additional metadata including 'cell_type', 'sex', and 'genotype'.

## Basic Usage:
### 1. Cross-species alignment
```
python ./src/cavebear_pytorch_cvae.py --input_h5ad ./data/example.h5ad --predict train
```
For hyperparameter tuning, we recommend tuning the learning_rate (--learining-rate {0.01 0.001 0.0001}).  Other potential hyperparameters to test are nlayers and embed_dim.  
> Defaults:  
&nbsp;&nbsp;--learning_rate 0.001  
&nbsp;&nbsp;--nlayers 3  
&nbsp;&nbsp;--embed_dim 25  

When predicting pseudotime across experimental systems (i.e. __in vitro_ and _in vivo_), include the discriminator argument and tune the discriminator weight.
> Arguments to include for discriminator:  
&nbsp;&nbsp;--dis dis   
&nbsp;&nbsp;--discriminator_weight {1.0 2.0 5.0 10.0} 


### 2. Cross-species model selection
Run 'get_best_params.py' with the path to the LISI_log.txt to select the model with optimal hyperparameter settings.
This creates a 'best_params.json' file in the same folder which can be used to set the arguments for pseudotime training and prediction.
```
python ./src/get_best_params.py --log /path/to/LISI_log.txt
```


### 3. Psuedotime training on reference cells and prediction on target cells
Use the 'best_params.json' file to get pseudotime prediction.
```
python ./src/cavebear_pytorch_cvae.py --input_h5ad ./data/example.h5ad --predict predict --train_species mouse --target_species zebrafish --time_label mouse_age
```
Where  `time_label` is the name of the column containing the training species time data (ie mouse age when cells were collected).

Output:
A .txt file with the model name and best hyperparameters for both trainings is created for the target cells by saving the adata.obs where the last column is called "pred_time" and contains the predicted pseudotime.
```
sampleID        age  batch  genotype  major_trajectory  origin  species  pred_time
GAP13.48.P9_G1	48.0000  0  XX      retinal neuron    Trapnell	zebrafish	  14.4456
GAP13.48.P9_C4	48.0000  0  XX      hatching gland    Trapnell	zebrafish	  16.2055
GAP13.48.P9_F8	48.0000  0  XX      hatching gland    Trapnell	zebrafish	  16.5047
GAP13.48.P9_H3	48.0000  0  XX      hatching gland    Trapnell	zebrafish	  14.3614
```


### 4. (OPTIONAL) Extract updated gene values after cross-species alignment
```
python /src/cavebear_pytorch_cvae.py --input_h5ad /data/example.h5ad --predict px_decoder 
```
Use best parameters from 'best_params.json' to ensure correct model is used.
Outputs are a species-agnostic gene expression numpy array, a metadata file (.tsv), and gene names (.txt).

# Cavebear: reference-guided pseudotime inference across species and biological contexts

<img width="1979" height="988" alt="Cavebear_overview_github" src="https://github.com/user-attachments/assets/7d81cbf9-1f39-45ce-8b3f-23335358b540" />

Overview: Cavebear is a reference-guided framework for predicting cellular pseudotime in undercharacterized datasets using time-series single-cell RNA-seq references from well-characterized systems.

Cavebear enables pseudotime inference across:

- species (e.g. mouse → human),
- experimental systems (e.g. _in vivo_ → _in vitro_),
- biological conditions (e.g. normal tissue → disease).

Method: Cavebear aligns species or experimental systems (step 1) and uses the learned cell embeddings of the reference cells to train the time predictor (step 2). The trained time predictor is subsequently applied to the learned cell embeddings of the target species or system from step 1 to predict the pseudotime of the target cells (step 3). 


## Installation:
Install through conda:
```
conda env create -f environment.yml
conda activate cavebear
```

## Input Data:
Cavebear requires a `.h5ad` file containing both reference and target datasets. Genes should be matched across datasets using ortholog mapping for cross-species applications or gene name matching for datasets from the same species.

Metadata should be stored in `adata.obs`.

### Required columns

| Column | Description |
|----------|------------|
| `species` | Labels identifying the reference and target datasets. These can correspond to different species, experimental systems, biological conditions, or datasets. |
| `batch` | Batch information. |
| `time` | Collection time of cells (e.g. developmental age, chronological age, or sampling time), only time labels in the reference data is used during training |

### Optional columns

- `cell_type` | this is not used in model training and may only be used in downstream evaluation/plotting
- `seed` | this can be used to to run model with various seeds (default is 101)


## Example Run:
```
bash ./src/run.sh --input {input.h5ad} --train_species {species1} --target_species {species2} [OPTIONAL ARGUMENTS]
```


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

When predicting pseudotime across experimental systems (e.g. _in vitro_ and _in vivo_), you may include the discriminator:
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
python ./src/cavebear_pytorch_cvae.py --input_h5ad ./data/example.h5ad --predict predict --train_species mouse --target_species zebrafish --time mouse_age
```
Where  `time` is the name of the column containing the training species time data (ie mouse age when cells were collected).

Output:
A .txt file with the model name and best hyperparameters for both trainings is created for the target cells by saving the adata.obs where the last column is called 'pred_time' and contains the predicted pseudotime.
```
sampleID        age  batch  genotype  major_trajectory  origin  species  pred_time
GAP13.48.P9_G1	48.0000  0  XX      retinal neuron    Trapnell	zebrafish	  14.4456
GAP13.48.P9_C4	48.0000  0  XX      hatching gland    Trapnell	zebrafish	  16.2055
GAP13.48.P9_F8	48.0000  0  XX      hatching gland    Trapnell	zebrafish	  16.5047
GAP13.48.P9_H3	48.0000  0  XX      hatching gland    Trapnell	zebrafish	  14.3614
```

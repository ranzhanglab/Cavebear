# Cross-species pseudotime prediction

### Overview: Cavebear utilizes a pytorch framework with two-steps to predict pseudotime of a sample using a well defined reference of a different species or experimental system (ie _in vivo_ vs _in vitro_).  

## Installation:
Install through conda:
```
conda env create -f environment.yml
conda activate
```

## Example Run:
```
cd src/
bash ./run.sh
```

## Basic Usage:
### 1. Cross-species alignment
```
python /src/cavebear_pytorch_cvae.py --input_h5ad {path/to/input_file.h5ad} --predict train
```
example input file: /data/example.h5ad
#### For hyperparameter tuning, the two parameters we recommend tuning are the learning_rate (--learining-rate {0.01 0.001 0.0001}) and discriminator (--dis dis & --discriminator_weight {1 2 5 10}). Other potential hyperparameters to test are nlayers and embed_dim.  
> Defaults:  
&nbsp;&nbsp;learning_rate = 0.001  
&nbsp;&nbsp;nlayers = 3  
&nbsp;&nbsp;embed_dim = 25  
&nbsp;&nbsp;dis=''  
&nbsp;&nbsp;discriminator_weight = 1.0  

The best model is selected by the highest LISI value which can be found in *_LISI.txt


### 2. Psuedotime training and prediction
```
python /src/cavebear_pytorch_cvae.py --input_h5ad /data/example.h5ad --predict predict --train_species mouse --target_species zebrafish --time_label mouse_age {optional: arguments for hyperparameters of best model if not default}
```
Where  `time_label` is the name of the column containing the training species time data (ie mouse age when cells were collected).  
If the best model from cross-species alignment has hyperparameters that differ from the defaults, you must include them as arguments.  


### 3. Extract gene probabilites after cross-species alignment
```
python /src/cavebear_pytorch_cvae.py --input_h5ad /data/example.h5ad --predict px_decoder {optional: arguments for hyperparameters of best model if not default}
```

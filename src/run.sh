#!/bin/bash

cur_dir=$(pwd)
nlayer=3
ndim=25

input_h5ad='../data/example.h5ad'
train_species=mouse
target_species=zebrafish

## 1. train the model
for lr in 0.01 0.001 0.0001; do
    python ${cur_dir}/cavebear_pytorch_cvae.py --input_h5ad ${input_h5ad} --predict train --learning_rate ${lr} --nlayer ${nlayer} --embed_dim ${ndim} --train_species ${train_species} --target_species ${target_species}
done


## 2. Select the best parameters -- outputs json file with best params (best_params.json)
python get_best_params.py --log /path/to/LISI_log.txt


## 3. Parse the json file to set best parameters and run the time predictor
BEST_PARAMS="/path/to/best_params.json"

LR=$(jq -r '.lr' $BEST_PARAMS)
N_LAYERS=$(jq -r '.n_layers' $BEST_PARAMS)
LATENT_DIM=$(jq -r '.latent_dim' $BEST_PARAMS)
DIS=$(jq -r '.dis' $BEST_PARAMS) # this is 0.0 if --dis was not used in the model, include '--dis' and '--discriminator_weight ${DIS}' arguments if used for selected best model

python ${cur_dir}/cavebear_pytorch_cvae.py --input_h5ad ${input_h5ad} --predict predict --learning_rate ${LR} --nlayer ${N_LAYERS} --embed_dim ${LATENT_DIM} --train_species ${train_species} --target_species ${target_species}

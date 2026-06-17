#!/bin/bash

# --- Arguments (input, train_species, and target_species are required) ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input | -i) input=$2; shift 2 ;;
        --train_species | -x) train_species=$2; shift 2 ;;
        --target_species | -y) target_species=$2; shift 2 ;;
        --nlayer | -l) nlayer=$2; shift 2 ;;
        --ndim | -d) ndim=$2; shift 2 ;;
        --extract | -e) extract=$2; shift 2 ;;
        --help   | -h)
            echo "Usage: $0 --input FILE --train_species STR --target_species STR [--nlayer INT --ndim INT --extract {true or false}]"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Defaults (applied only if flag was not provided) ---
nlayer=${nlayer:-3}
ndim=${ndim:-25}
extract=${ndim:-false}

# --- Validate arguments ---
error=false

if [[ -z "$input" ]]; then
    echo "Error: --input is required"
    error=true
fi

if [[ -z "$train_species" ]]; then
    echo "Error: --training_species is required"
    error=true
fi

if [[ -z "$target_species" ]]; then
    echo "Error: --target_species is required"
    error=true
fi

if [[ "$error" == true ]]; then
    echo "Usage: $0 --input FILE --train_species STR --target_species STR [--learning_rate FLOAT --nlayer INT --ndim INT]"
    exit 1
fi


## --- Set directory and get sample basename ----------------------------------------------------------
cur_dir=$(pwd)/src
input_base=$(basename "$input")
name="${input_base%.*}"



## --- 1. train the model ----------------------------------------------------------
for lr in 0.01 0.001 0.0001; do
    python ${cur_dir}/cavebear_pytorch_cvae.py --input_h5ad ${input} --predict train --learning_rate ${lr} --nlayer ${nlayer} --embed_dim ${ndim} --train_species ${train_species} --target_species ${target_species}
done


## --- 2. Select the best parameters -- outputs json file with best params (best_params.json) ----------------------------------------------------------
python get_best_params.py --log ../${cur_dir}/results/${name}/LISI_log.txt


## --- 3. Parse the json file to set best parameters and run the time predictor ----------------------------------------------------------
BEST_PARAMS="../${cur_dir}/results/${name}/best_params.json"

LR=$(jq -r '.lr' $BEST_PARAMS)
N_LAYERS=$(jq -r '.n_layers' $BEST_PARAMS)
LATENT_DIM=$(jq -r '.latent_dim' $BEST_PARAMS)
DIS=$(jq -r '.dis' $BEST_PARAMS) # this is 0.0 if --dis was not used in the model, include '--dis' and '--discriminator_weight ${DIS}' arguments if used for selected best model

python ${cur_dir}/cavebear_pytorch_cvae.py --input_h5ad ${input} --predict predict --learning_rate ${LR} --nlayer ${N_LAYERS} --embed_dim ${LATENT_DIM} --train_species ${train_species} --target_species ${target_species}


## --- 4. (Optional) Extract gene expression values after cross-species alignment (step 1) ----------------------------------------------------------
if [[ "$extract" == true ]]; then
    python ${cur_dir}/cavebear_pytorch_cvae.py --input_h5ad ${input} --predict px_decoder --learning_rate ${LR} --nlayer ${N_LAYERS} --embed_dim ${LATENT_DIM} --train_species ${train_species} --target_species ${target_species}
fi
#!/bin/bash

# --- Arguments (input, train_species, and target_species are required) ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input | -i) input=$2; shift 2 ;;
        --train_species | -x) train_species=$2; shift 2 ;;
        --target_species | -y) target_species=$2; shift 2 ;;
        --cell_type | -c) cell_type=$2; shift 2 ;;
        --time | -t) time=$2; shift 2 ;;
        --use_dis | -u) use_dis=$2; shift 2 ;;
        --seed | -s) seed=$2; shift 2 ;;
        --help   | -h)
            echo "Usage: $0 --input FILE --train_species STR --target_species STR --cell_type STR --time {STR; default: time} --use_dis {true or false; default: false} --seed {INT; default: 101}]"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Defaults (applied only if flag was not provided) ---
cell_type=${cell_type:-''}
time=${time:-'time'}
use_dis=${use_dis:-false}
seed=${seed:-101}

# --- Validate arguments ---
error=false

if [[ -z "$input" ]]; then
    echo "Error: --input is required"
    error=true
fi

if [[ ! -f "$input" ]]; then
    echo "Error: input file not found: $input"
    error=true
fi

if [[ -z "$train_species" ]]; then
    echo "Error: --train_species is required"
    error=true
fi

if [[ -z "$target_species" ]]; then
    echo "Error: --target_species is required"
    error=true
fi

if [[ "$error" == true ]]; then
    echo "Usage: $0 --input FILE --train_species STR --target_species STR --cell_type STR --time {STR; default: time} --use_dis {true or false; default: false} --seed {INT; default: 101}]"
    exit 1
fi


## --- Set directory and get sample basename ----------------------------------------------------------
cur_dir=$(pwd)
script_dir=$(pwd)/src
input_base=$(basename "$input")
name="${input_base%.*}"



## --- 1. train the model ----------------------------------------------------------
cd ${script_dir}

# set seed argument if not default
if [[ "$seed" != "101" ]]; then
    SEED_ARGS="--seed ${seed}"
else
    SEED_ARGS=""
fi

if [[ "$use_dis" == "false" ]]; then
    for lr in 0.01 0.001 0.0001; do
        python ${script_dir}/cavebear_pytorch_cvae.py \
        --input_h5ad ${input} \
        --predict train \
        --learning_rate ${lr} \
        --train_species ${train_species} \
        --target_species ${target_species} \
        ${SEED_ARGS}
    done
else
    for lr in 0.01 0.001 0.0001; do
        for dis_weight in 0 1 2 5 10 20; do
            python ${script_dir}/cavebear_pytorch_cvae.py \
            --input_h5ad ${input} \
            --predict train \
            --learning_rate ${lr} \
            --train_species ${train_species} \
            --target_species ${target_species} \
            --dis dis \
            --discriminator_weight ${dis_weight} \
            ${SEED_ARGS}
        done
    done
fi

## --- 2. Select the best parameters -- outputs json file with best params (best_params.json) ----------------------------------------------------------
python ${script_dir}/get_best_params.py --log ${cur_dir}/results/${name}/LISI_log.txt


## --- 3. Parse the json file to set best parameters and run the time predictor ----------------------------------------------------------
BEST_PARAMS=${cur_dir}/results/${name}/best_params.json

# check that best_params.json exists
if [[ ! -f "$BEST_PARAMS" ]]; then
    echo "Error: best params file not found: $BEST_PARAMS"
    exit 1
fi

LR=$(jq -r '.lr' $BEST_PARAMS)
N_LAYERS=$(jq -r '.n_layers' $BEST_PARAMS)
LATENT_DIM=$(jq -r '.latent_dim' $BEST_PARAMS)
DIS=$(jq -r '.dis' $BEST_PARAMS) 
SEED=$(jq -r '.seed' $BEST_PARAMS)

# set DIS_ARGS and PARAM_STRING (used in step 4)
if [[ "$use_dis" == "true" ]]; then
    DIS_ARGS="--dis dis --discriminator_weight ${DIS}"
    PARAM_STRING="${LR}_${N_LAYERS}_${LATENT_DIM}_dis${DIS}"
else
    DIS_ARGS=""
    PARAM_STRING="${LR}_${N_LAYERS}_${LATENT_DIM}"
fi

if [[ "$SEED" != "101" ]]; then
    SEED_ARGS="--seed ${SEED}"
    PARAM_STRING="${PARAM_STRING}_seed${SEED}"
else
    SEED_ARGS=""
fi

python ${script_dir}/cavebear_pytorch_cvae.py \
    --input_h5ad ${input} \
    --predict time \
    --learning_rate ${LR} \
    --train_species ${train_species} \
    --target_species ${target_species} \
    --time ${time} \
    ${DIS_ARGS} \
    ${SEED_ARGS}


## -- 4. Evaluate the time prediction by pairwise accuracy -- option to split by cell_type
model_input="${cur_dir}/results/${name}/${PARAM_STRING}/cvae_pytorch_disc_best_model_${name}_${PARAM_STRING}.pth"

# check that model_input exists
if [[ ! -f "$model_input" ]]; then
    echo "Error: model_input file not found: $model_input"
    exit 1
fi


if [ -n "$cell_type" ]; then
    Rscript ${script_dir}/eval_time_pred.R \
        --input ${model_input} \
        --target_species ${target_species} \
        --time ${time} \
        --cell_type ${cell_type}
else
        Rscript ${script_dir}/eval_time_pred.R \
        --input ${model_input} \
        --target_species ${target_species} \
        --time ${time}
fi


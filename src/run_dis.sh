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
        --cell_type | -e) cell_type=$2; shift 2 ;;
        --time_label | -e) time_label=$2; shift 2 ;;
        --help   | -h)
            echo "Usage: $0 --input FILE --train_species STR --target_species STR [--nlayer INT --ndim INT --extract {true or false}]"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Defaults (applied only if flag was not provided) ---
nlayer=${nlayer:-3}
ndim=${ndim:-25}
extract=${extract:-false}
cell_type=${cell_type:-''}
time_label=${time_label:-'time'}

# --- Validate arguments ---
error=false

if [[ -z "$input" ]]; then
    echo "Error: --input is required"
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
    echo "Usage: $0 --input FILE --train_species STR --target_species STR [--learning_rate FLOAT --nlayer INT --ndim INT]"
    exit 1
fi


## --- Set directory and get sample basename ----------------------------------------------------------
cur_dir=$(pwd)
script_dir=$(pwd)/src
input_base=$(basename "$input")
name="${input_base%.*}"



## --- 1. train the model ----------------------------------------------------------
cd ${script_dir}
for lr in 0.01 0.001 0.0001; do
    for dis_weight in 0 1 2 5 10; do
        python ${script_dir}/cavebear_pytorch_cvae.py \
        --input_h5ad ${input} \
        --predict train \
        --learning_rate ${lr} \
        --nlayer ${nlayer} \
        --embed_dim ${ndim} \
        --train_species ${train_species} \
        --target_species ${target_species} \
        --dis dis \
        --discriminator_weight ${dis_weight}
    done
done


## --- 2. Select the best parameters -- outputs json file with best params (best_params.json) ----------------------------------------------------------
python ${script_dir}/get_best_params.py --log ${cur_dir}/results/${name}/LISI_log.txt


## --- 3. Parse the json file to set best parameters and run the time predictor ----------------------------------------------------------
BEST_PARAMS=${cur_dir}/results/${name}/best_params.json

LR=$(jq -r '.lr' $BEST_PARAMS)
N_LAYERS=$(jq -r '.n_layers' $BEST_PARAMS)
LATENT_DIM=$(jq -r '.latent_dim' $BEST_PARAMS)
DIS=$(jq -r '.dis' $BEST_PARAMS) 

if (( $(echo "$DIS != 0.0" | bc -l) )); then
    DIS_ARGS="--dis dis --discriminator_weight ${DIS}"
    PARAM_STRING="${LR}_${N_LAYERS}_${LATENT_DIM}_dis${DIS}"
else
    DIS_ARGS="--dis dis --discriminator_weight ${DIS}"
    PARAM_STRING="${LR}_${N_LAYERS}_${LATENT_DIM}_dis"
fi

python ${script_dir}/cavebear_pytorch_cvae.py \
    --input_h5ad ${input} \
    --predict time \
    --learning_rate ${LR} \
    --nlayer ${N_LAYERS} \
    --embed_dim ${LATENT_DIM} \
    --train_species ${train_species} \
    --target_species ${target_species} \
    --time_label ${time_label} \
    ${DIS_ARGS}


## -- 4. Evaluate the time prediction by pairwise accuracy -- option to split by cell_type
if [ "$target_species" == "human" ]; then
    PARAM_STRING="${PARAM_STRING}_human"
fi

model_input="${cur_dir}/results/${name}/${PARAM_STRING}/cvae_pytorch_disc_best_model_${name}_${PARAM_STRING}.pth"

if [ -n "$cell_type" ]; then
    Rscript ${script_dir}/eval_time_pred.R \
        --input ${model_input} \
        --target_species ${target_species} \
        --time ${time_label} \
        --cell_type ${cell_type}
else
        Rscript ${script_dir}/eval_time_pred.R \
        --input ${model_input} \
        --target_species ${target_species} \
        --time ${time_label}
fi


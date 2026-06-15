#!/bin/bash

cur_dir=$(pwd)
ndim = 25

## train the model
for lr in 0.01 0.001 0.0001; do
    for nlayer in 2 3 4; do
        python ${cur_dir}/cavebear_pytorch_cvae.py --learning_rate ${lr} --nlayer ${nlayer} --embed_dim ${ndim} --predict train --train_species mouse --target_species zebrafish
    done
done



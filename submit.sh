#!/bin/bash

for lambda in $(seq 0 0.5 9)
do
    sed -i "s/lambda=.*/lambda=$lambda/" input.nml
    make run_dyn
    make run_dyn_g
done

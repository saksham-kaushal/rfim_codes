#!/bin/bash

cp -r ../rfim_correlated/data/Phi-0-1800.csv ./data_input
cp -r ../rfim_correlated/data/Si-0-{16..20}00.csv ./data_input
cp beta.txt ./data_input

rename 'Si-0-' 'Si-0_0-' ./data_input/Si*.csv
rename 'Phi-0-' 'Phi-0-0-' ./data_input/Phi*.csv

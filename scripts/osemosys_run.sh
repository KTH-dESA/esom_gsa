#!/usr/bin/env bash
MODELRUN=$1
RESULTS="processed_data\/$MODELRUN\/SelectedResults.csv"
mkdir processed_data/$MODELRUN
cat model/osemosys.txt > processed_data/$MODELRUN/osemosys.txt
sed -i '' "s/FILEPATH/$RESULTS/g" processed_data/$MODELRUN/osemosys.txt

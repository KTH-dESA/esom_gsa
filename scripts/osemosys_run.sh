#!/usr/bin/env bash
MODELRUN=$1
RESULTS="processed_data\/$MODELRUN\/SelectedResults.csv"
mkdir processed_data/$MODELRUN
cat model/osemosys_short.txt > processed_data/$MODELRUN/osemosys_short.txt
sed -i '' "s/FILEPATH/$RESULTS/g" processed_data/$MODELRUN/osemosys_short.txt

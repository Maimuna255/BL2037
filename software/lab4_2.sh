#!/bin/bash

# Definiera sökvägar i variabler (valfritt men snyggt)
INPUT="./data/lab3/gene_expression.csv"
OUTPUT="./results/lab4/hist.png"

# Kör python-skriptet med argumenten
python3 lab4_1.py $INPUT $OUTPUT

#!/bin/bash

# Skapa resultatmappen om den inte finns
mkdir -p ./results/lab4

# 2. Kör första filen
# Vi pekar på software/lab4_1.py och data/lab3/...
python3 software/lab4_1.py ./data/lab3/gene_expression.csv ./results/lab4/hist1.png &

# 3. Kör andra filen (om gene_expression2.csv också ligger i data/lab3)
python3 software/lab4_1.py ./data/lab3/gene_expression2.csv ./results/lab4/hist2.png &

# Vänta på att båda blir klara
wait

echo "Båda analyserna är klara! Bilderna finns i ./results/lab4/"

#!/bin/bash

# Skapa resultatmappen om den inte finns
mkdir -p ./results/lab4

# Kör första filen i bakgrunden (&)
python3 lab4_1.py ./examples/gene_expression.csv ./results/lab4/hist1.png &

# Kör andra filen i bakgrunden (&)
python3 lab4_1.py ./examples/gene_expression2.csv ./results/lab4/hist2.png &

# Vänta på att båda processerna ska bli klara innan skriptet avslutas
wait

echo "Båda analyserna är klara! Bilderna finns i ./results/lab4/"

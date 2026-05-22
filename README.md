# Protein Pipeline Project

Det här projektet innehåller en automatiserad bioinformatik-pipeline skriven i Python för att analysera och jämföra proteinstrukturer från PDB och AlphaFold.

## Innehåll
* `pipeline.py` - Hämtar experimentella data från PDB, extraherar sekvenser och kör Clustal Omega.
* `alphafold_pipeline.py` - Analyserar motsvarande sekvenser utifrån ett AlphaFold-perspektiv.
* `inputs_finalproject.csv` - Inputdata med de 25 målproteinerna.

## Hur man kör koden
1. Kör PDB-analysen: `python3 finalprojectscript/pipeline.py`
2. Kör AlphaFold-analysen: `python3 finalprojectscript/alphafold_pipeline.py`

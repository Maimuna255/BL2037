import pandas as pd
import Bio
from Bio.PDB import MMCIFParser, PDBList
import os
import subprocess

print("--- STARTAR STEG 1: Laddar ner CIF-filer ---")

csv_fil = 'inputs_finalproject.csv'

if not os.path.exists(csv_fil):
    print(f"FEL: Hittade inte filen {csv_fil} i den här mappen!")
else:
    df = pd.read_csv(csv_fil)
    print(f"Lyckades läsa CSV! Hittade {len(df)} rader.")

    if not os.path.exists('finalprojectdata'):
        os.makedirs('finalprojectdata')

    pdb_kolumn = df.columns[0] 
    pdb_ids = df[pdb_kolumn].dropna().unique()

    pdbl = PDBList()
    for pid in pdb_ids:
        pid = str(pid).strip()
        pdb_id_clean = pid.split('_')[0]
        pdbl.retrieve_pdb_file(pdb_id_clean, pdir='finalprojectdata', file_format='mmCif')

    print("\n--- STEG 1 KLART: Alla CIF-filer har laddats ner! ---")

    print("\n--- STARTAR STEG 2: Extraherar sekvenser från CIF till FASTA ---")
    from Bio.SeqUtils import seq1

    fasta_filnamn = "finalprojectdata/sequences.fasta"

    with open(fasta_filnamn, "w") as fasta_fil:
        parser = MMCIFParser(QUIET=True)
        
        for _, row in df.iterrows():
            full_id = str(row[pdb_kolumn]).strip()
            parts = full_id.split('_')
            pdb_id = parts[0].lower()
            chain_id = parts[1] if len(parts) > 1 else 'A'
            
            lokal_fil = f"finalprojectdata/{pdb_id}.cif"
            if not os.path.exists(lokal_fil):
                continue
                
            try:
                structure = parser.get_structure(pdb_id, lokal_fil)
                model = structure[0]
                
                # Speciallösning för 4E5A där kedjan heter X istället för A
                if pdb_id == "4e5a" and "X" in model:
                    chain_id = "X"
                
                if chain_id in model:
                    chain = model[chain_id]
                    sekvens_bokstäver = ""
                    for residue in chain:
                        if residue.get_id()[0] == " ": 
                            sekvens_bokstäver += seq1(residue.get_resname())
                    
                    sekvens_bokstäver = sekvens_bokstäver.replace('X', '').replace('?', '')
                    
                    if len(sekvens_bokstäver) > 0:
                        fasta_fil.write(f">{full_id}\n{sekvens_bokstäver}\n")
                        print(f"Extraherade kedja {chain_id} från {pdb_id.upper()} ({len(sekvens_bokstäver)} aa)")
                else:
                    print(f"Hoppar över {full_id} (kedja {chain_id} saknas)")
                    
            except Exception as e:
                print(f"Kunde inte processa {pdb_id.upper()}: {e}")

    print(f"\n--- STEG 2 KLART: FASTA-fil skapad! ---")

    print("\n--- STARTAR STEG 3: Kör Multiple Sequence Alignment med Clustal Omega ---")
    
    output_alignment = "finalprojectdata/aligned_sequences.fasta"
    
    # Kommandot för att köra clustalo via terminalen
    clustal_cmd = [
        "clustalo",
        "-i", fasta_filnamn,
        "-o", output_alignment,
        "--outfmt=fa",
        "--force"
    ]
    
    try:
        print("Kör Clustal Omega...")
        subprocess.run(clustal_cmd, check=True)
        print(f"ALIGNMENT KLART! Resultatet har sparats i: {output_alignment}")
    except FileNotFoundError:
        print("FEL: Kommandot 'clustalo' hittades inte. Har du installerat det via Conda/Brew?")
    except subprocess.CalledProcessError as e:
        print(f"FEL vid körning av Clustal Omega: {e}")

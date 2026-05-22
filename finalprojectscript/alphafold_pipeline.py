import pandas as pd
import Bio
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import os
import subprocess
import matplotlib.pyplot as plt

print("--- STARTAR STEG 6 (BETYG C): ALPHAFOLD / PDB JÄMFÖRELSE (LOKAL VERSION) ---")

csv_fil = 'inputs_finalproject.csv'
pdb_mapp = 'finalprojectdata'
data_mapp = 'alphafold_data'

if not os.path.exists(data_mapp):
    os.makedirs(data_mapp)

if not os.path.exists(csv_fil):
    print(f"FEL: Hittade inte filen {csv_fil}!")
    exit()

df = pd.read_csv(csv_fil)
pdb_kolumn = df.columns[0]

print("\nLäser sekvenser lokalt från dina nedladdade PDB CIF-filer...")
fasta_filnamn = f"{data_mapp}/af_sequences.fasta"
parser = MMCIFParser(QUIET=True)

har_data = False
with open(fasta_filnamn, "w") as fasta_fil:
    for _, row in df.iterrows():
        full_id = str(row[pdb_kolumn]).strip()
        parts = full_id.split('_')
        pdb_id = parts[0].lower()
        
        # Speciallösning från tidigare steg för att läsa 4E5A
        chain_id = 'X' if pdb_id == '4e5a' else (parts[1] if len(parts) > 1 else 'A')
        
        lokal_fil = f"{pdb_mapp}/{pdb_id}.cif"
        if not os.path.exists(lokal_fil):
            continue  # Hoppa över de som saknas helt (t.ex. 4v88, 6z6l)
            
        try:
            structure = parser.get_structure(pdb_id, lokal_fil)
            model = structure[0]
            
            if chain_id in model:
                chain = model[chain_id]
                sekvens = "".join([seq1(r.get_resname()) for r in chain if r.get_id()[0] == " "])
                sekvens = sekvens.replace('X', '').replace('?', '')
                
                if len(sekvens) > 0:
                    fasta_fil.write(f">{full_id}\n{sekvens}\n")
                    print(f"  -> Extraherade {len(sekvens)} aa för {full_id}")
                    har_data = True
        except Exception as e:
            print(f"  -> Fel vid läsning av {full_id}: {e}")

if not har_data:
    print("\n[FEL]: Kunde inte läsa lokala filer. Körde du Steg 1 först?")
    exit()

print("\nKör Multiple Sequence Alignment med Clustal Omega...")
aligned_fasta = f"{data_mapp}/af_aligned.fasta"
clustal_cmd = ["clustalo", "-i", fasta_filnamn, "-o", aligned_fasta, "--outfmt=fa", "--force"]

try:
    subprocess.run(clustal_cmd, check=True)
    print(f"  -> Alignment klart!")
except Exception as e:
    print(f"  -> FEL vid alignment: {e}")
    exit()

print("\nBygger fylogenetiskt träd...")
try:
    alignment = AlignIO.read(aligned_fasta, "fasta")
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)

    for clade in tree.find_clades():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = ""

    fig = plt.figure(figsize=(10, 8), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.title("Fylogenetiskt träd (AlphaFold-sekvensanalys)", fontsize=12, fontweight='bold')
    plt.xlabel("Evolutionärt avstånd")
    plt.ylabel("Proteiner")
    plt.tight_layout()
    
    output_bild = f"{data_mapp}/alphafold_tree.png"
    plt.savefig(output_bild)
    print(f"\n--- ALLT KLART! AlphaFold-träd sparat i: {output_bild} ---")
except Exception as e:
    print(f"FEL vid trädbygge: {e}")

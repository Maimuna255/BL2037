import subprocess
from Bio import PDB, SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import io

def get_sequence_from_pdb(pdb_id, file_path):
    """Extraherar sekvensen för kedja R från en PDB-fil."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, file_path)
    
    # Tre-bokstavskod till en-bokstavskod
    ppb = PDB.PPBuilder()
    seq = ""
    for pp in ppb.build_peptides(structure[0]['R']):
        seq += str(pp.get_sequence())
    return seq

def main():
    pdb_files = {
        "8JIQ_R": "8JIQ_R.pdb",
        "8E3Y_R": "8E3Y_R.pdb",
        "6WI9_R": "6WI9_R.pdb",
        "6X18_R": "6X18_R.pdb",
        "7VQX_R": "7VQX_R.pdb",
        "7YON_R": "7YON_R.pdb"
    }

    # 1. Skapa FASTA-fil
    fasta_records = []
    for name, path in pdb_files.items():
        sequence = get_sequence_from_pdb(name, path)
        record = SeqRecord(Seq(sequence), id=name, description="")
        fasta_records.append(record)
    
    with open("receptors.fasta", "w") as f:
        SeqIO.write(fasta_records, f, "fasta")
    print("✓ receptors.fasta har skapats.")

    # 2. Kör Clustal Omega via subprocess
    input_file = "receptors.fasta"
    fasta_output = "aligned_receptors.fasta"
    tree_output = "tree_receptors.nwk"
    
    clustal_command = [
        "clustalo", 
        "-i", input_file, 
        "-o", fasta_output, 
        "--outfmt=fasta", 
        "--guidetree-out", tree_output, 
        "--force"
    ]

    print("Kör Clustal Omega...")
    subprocess.run(clustal_command, check=True)
    print("✓ Sekvensjustering klar.")

    # 3. Läs och visa trädet i terminalen
    print("\n--- Fylogenetiskt träd (Sequence Similarity) ---")
    tree = Phylo.read(tree_output, "newick")
    Phylo.draw_ascii(tree)

if __name__ == "__main__":
    main()

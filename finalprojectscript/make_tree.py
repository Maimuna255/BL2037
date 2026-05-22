import os
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

print("--- STARTAR TRÄDBYGGANDE ---")

alignment_fil = "finalprojectdata/aligned_sequences.fasta"
output_bild = "finalprojectdata/phylogenetic_tree.png"

if not os.path.exists(alignment_fil):
    print(f"FEL: Hittade inte alignment-filen {alignment_fil}!")
else:
    # 1. Läs in din Multiple Sequence Alignment
    print("Läser in alignment...")
    alignment = AlignIO.read(alignment_fil, "fasta")

    # 2. Beräkna avståndsmatris (Distance Matrix) baserat på aminosyralikhet
    print("Räknar ut evolutionära avstånd (Identity-modell)...")
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    # 3. Bygg trädet med Neighbor Joining-metoden (en standardmetod för proteiner)
    print("Konstruerar fylogenetiskt träd (Neighbor Joining)...")
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)

    # 4. Snygga till trädet (ta bort interna nodnamn som kan se stökiga ut)
    for clade in tree.find_clades():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = ""

    # 5. Rita och spara trädet som en bild
    print("Genererar och sparar bild...")
    fig = plt.figure(figsize=(10, 8), dpi=300) # Gör bilden stor och högupplöst
    ax = fig.add_subplot(1, 1, 1)
    
    # Rita trädet via Biopythons inbyggda funktion
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Justera designen lite
    plt.title("Fylogenetiskt träd baserat på proteinsekvenser", fontsize=14, fontweight='bold')
    plt.xlabel("Evolutionärt avstånd", fontsize=11)
    plt.ylabel("Proteiner", fontsize=11)
    plt.tight_layout()
    
    # Spara filen
    plt.savefig(output_bild)
    print(f"\nKLART! Ditt fylogenetiska träd har sparats som en bild i: {output_bild}")

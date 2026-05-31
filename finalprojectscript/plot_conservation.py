import os
import matplotlib.pyplot as plt
from Bio import AlignIO
from collections import Counter

print("--- STARTAR GENERERING AV KONSERVERINGSPLOT ---")

alignment_fil = "finalprojectdata/aligned_sequences.fasta"
output_bild = "finalprojectdata/conservation_plot.png"

if not os.path.exists(alignment_fil):
    print(f"FEL: Hittade inte alignment-filen: {alignment_fil}")
else:
    # 1. Läs in din Multiple Sequence Alignment (MSA)
    alignment = AlignIO.read(alignment_fil, "fasta")
    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)
    
    conservation_scores = []
    
    # 2. Räkna ut likheten för varje enskild aminosyra-position
    for i in range(alignment_length):
        column = alignment[:, i]
        counts = Counter(column)
        
        # Hitta den absolut vanligaste aminosyran på denna specifika position
        most_common_count = counts.most_common(1)[0][1]
        
        # Räkna ut procentuell konservering (mellan 0.0 och 1.0)
        score = most_common_count / num_sequences
        conservation_scores.append(score)
        
    # 3. Skapa och rita linjegrafen
    plt.figure(figsize=(12, 5), dpi=300)
    plt.plot(conservation_scores, color="teal", linewidth=1.2, label="Konserveringsgrad")
    
    # Snygga till designen för rapporten
    plt.title("Sequence Conservation Score per Alignment Position", fontsize=14, fontweight='bold')
    plt.xlabel("Alignment Position (Aminosyra-index)", fontsize=11)
    plt.ylabel("Conservation Score (0.0 = Variabel, 1.0 = Helt konserverad)", fontsize=11)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.ylim(-0.05, 1.05)
    plt.tight_layout()
    
    # 4. Spara bilden
    plt.savefig(output_bild)
    print(f"\nKLART! Din konserveringsgraf har sparats i: {output_bild}")

import subprocess
import re
import csv

def parse_tmalign_output(output_text):
    # Regex för att hitta TM-score och RMSD
    tm_pattern = r"TM-score\s*=\s*(\d\.\d+)"
    rmsd_pattern = r"RMSD\s*=\s*(\d\.\d+)"
    
    tm_match = re.search(tm_pattern, output_text)
    rmsd_match = re.search(rmsd_pattern, output_text)
    
    tm_score = tm_match.group(1) if tm_match else "N/A"
    rmsd = rmsd_match.group(1) if rmsd_match else "N/A"
    
    return tm_score, rmsd

def main():
    reference = "/Users/monirulislam/bioprojekt/structures_R/8JIQ_R.pdb"
    targets = ["/Users/monirulislam/bioprojekt/structures_R/8E3Y_R.pdb","/Users/monirulislam/bioprojekt/structures_R/6WI9_R.pdb","/Users/monirulislam/bioprojekt/structures_R/6X18_R.pdb","/Users/monirulislam/bioprojekt/structures_R/7VQX_R.pdb","/Users/monirulislam/bioprojekt/structures_R/7YON_R.pdb"]
    results = []

    print(f"{'Target':<15} | {'TM-score':<10} | {'RMSD':<10}")
    print("-" * 40)

    for target_file in targets:
        output_prefix = f"align_{target_file.split('.')[0]}"
        # Kommandot som körs: TMalign ref target -o prefix
        command = ["TMalign", reference, target_file, "-o", output_prefix]

        process = subprocess.run(command, capture_output=True, text=True)
        tm_score, rmsd = parse_tmalign_output(process.stdout)
        
        results.append({"PDB_ID": target_file, "TM-score": tm_score, "RMSD": rmsd})
        print(f"{target_file:<15} | {tm_score:<10} | {rmsd:<10}")

    # Spara till CSV
    with open("tmalign_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["PDB_ID", "TM-score", "RMSD"])
        writer.writeheader()
        writer.writerows(results)

if __name__ == "__main__":
    main()

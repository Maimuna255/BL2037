import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def main():
    # 1. Skapa en parser-instans
    parser = argparse.ArgumentParser(description="Skapar ett histogram från genuttrycksdata.")

    # 2. Lägg till argumenten
    parser.add_argument("input_csv", help="Sökväg till CSV-filen")
    parser.add_argument("output_png", help="Sökväg till PNG-filen")

    # 3. Läs in argumenten
    args = parser.parse_args()

    # 4. Använd argumenten för att läsa och spara
    df = pd.read_csv(args.input_csv)
    
    plt.figure()
    # Här väljer vi kolumn 2 (index 1) för histogrammet
    plt.hist(df.iloc[:, 1], bins=20)
    plt.title("Genuttryck")

    # Skapa mappen om den inte finns
    output_dir = os.path.dirname(args.output_png)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.savefig(args.output_png)
    print(f"Filen har sparats som {args.output_png}")

if __name__ == "__main__":
    main()

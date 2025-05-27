import argparse
import re

def parse_fasta(input_file, threshold):
    above_threshold = []
    below_threshold = []
    
    tpm_pattern = re.compile(r"\btpm:([0-9.]+)\b")  # Regex to match "tpm:<value>"
    
    with open(input_file, 'r') as f:
        header = ""
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # New header line
                if header:  # Previous sequence exists, add it to the correct list
                    match = tpm_pattern.search(header)
                    if match:
                        tpm_value = float(match.group(1))
                        if tpm_value >= threshold:
                            above_threshold.append(f"{header}\n{sequence}")
                        else:
                            below_threshold.append(f"{header}\n{sequence}")
                
                header = line  # Start new header
                sequence = ""  # Reset sequence
            else:
                sequence += line  # Append sequence data

        # Add the last sequence
        if header:
            match = tpm_pattern.search(header)
            if match:
                tpm_value = float(match.group(1))
                if tpm_value >= threshold:
                    above_threshold.append(f"{header}\n{sequence}")
                else:
                    below_threshold.append(f"{header}\n{sequence}")

    return above_threshold, below_threshold

def write_fasta(output_file, fasta_entries):
    with open(output_file, 'w') as f:
        for entry in fasta_entries:
            f.write(entry + "\n")

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Separate sequences based on TPM threshold.")
    parser.add_argument("input_file", help="Input FASTA file.")
    parser.add_argument("threshold", type=float, help="TPM threshold.")
    
    args = parser.parse_args()

    # Parse the input file and separate sequences
    above_threshold, below_threshold = parse_fasta(args.input_file, args.threshold)
    
    # Write the separated sequences to the respective files
    write_fasta('aboveThreshold.fa', above_threshold)
    write_fasta('belowThreshold.fa', below_threshold)

if __name__ == "__main__":
    main()

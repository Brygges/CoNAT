import sys

def filter_fasta(input_file, min_length, output_file):
    total_sequences = 0
    filtered_count = 0

    with open(input_file, 'r') as file, open(output_file, 'w') as out_file:
        sequence_id = None
        sequence = ""

        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if sequence_id and len(sequence) >= min_length:
                    out_file.write(sequence_id + "\n")
                    out_file.write(sequence + "\n")
                    filtered_count += 1
                total_sequences += 1
                sequence_id = line
                sequence = ""
            else:  # Sequence line
                sequence += line

        # Check the last sequence after exiting the loop
        if sequence_id and len(sequence) >= min_length:
            out_file.write(sequence_id + "\n")
            out_file.write(sequence + "\n")
            filtered_count += 1

    return total_sequences, filtered_count

def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_fasta.py <input_fasta_file> <min_length> <output_fasta_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    try:
        min_length = int(sys.argv[2])
    except ValueError:
        print("Error: min_length must be an integer.")
        sys.exit(1)
    
    output_file = sys.argv[3]
    
    total_sequences, filtered_count = filter_fasta(input_file, min_length, output_file)
    print(f"Total sequences in input file: {total_sequences}")
    print(f"Filtered sequences written to {output_file}: {filtered_count}")

if __name__ == "__main__":
    main()

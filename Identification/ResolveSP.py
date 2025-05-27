import argparse

def process_3line(input_file, fasta_file, cleaned_file):
    # Set to track cleaned sequences we've already processed
    seen_cleaned_sequences = set()

    with open(input_file, 'r') as infile, \
         open(fasta_file, 'w') as fasta_out, \
         open(cleaned_file, 'w') as cleaned_out:

        # Initialize variables for reading the file
        lines = infile.readlines()
        i = 0
        while i < len(lines):
            header = lines[i].strip() # Read the header line
            sequence = lines[i+1].strip() # Read the sequence line
            annotation = lines[i+2].strip() # Read the annotation line
            i += 3

            # Ensure that the sequence and annotation are the same length
            if len(sequence) != len(annotation):
                print(f"Warning: Skipping entry with mismatched sequence and annotation length: {header}")
                continue

            # Only process entries with SP in the header
            if '| SP' in header:
                cleaned_sequence = sequence.rstrip('*')
                if not cleaned_sequence:
                    print(f"Warning: Skipping empty sequence for entry: {header}")
                    continue

                # Count the number of 'S' in the annotation to determine how much of the sequence to discard
                signal_peptide_length = annotation.count('S')

                # Remove the signal peptide (S's) from the sequence
                cleaned_sequence = cleaned_sequence[signal_peptide_length:]

                # If the cleaned sequence is empty after removing the signal peptide, skip it
                if not cleaned_sequence:
                    print(f"Warning: Skipping cleaned sequence for entry with no remaining sequence after signal peptide removal: {header}")
                    continue
                if cleaned_sequence in seen_cleaned_sequences:
                    continue
                seen_cleaned_sequences.add(cleaned_sequence)

                # Save the sequence to the .fasta file (with signal peptide)
                fasta_out.write(f"{header}\n{sequence}\n")

                # Save the cleaned sequence (signal peptide removed) to the cleaned file
                cleaned_out.write(f"{header}\n{cleaned_sequence}\n")


def main():
    # Set up the command line argument parser
    parser = argparse.ArgumentParser(description="Process a 3-line sequence file, remove signal peptides, and save to two separate FASTA files.")
    
    # Define the expected arguments
    parser.add_argument("input_file", help="The input 3-line format file")
    parser.add_argument("fasta_file", help="The output FASTA file with signal peptide sequences")
    parser.add_argument("cleaned_file", help="The output FASTA file with cleaned sequences (signal peptide removed)")

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the process function with the arguments
    process_3line(args.input_file, args.fasta_file, args.cleaned_file)


if __name__ == "__main__":
    main()

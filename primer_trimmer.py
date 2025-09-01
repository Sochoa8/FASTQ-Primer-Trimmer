#!/usr/bin/env python3

import sys
import re

try:
    import tkinter as tk
    from tkinter import filedialog
except ImportError:
    print("Error: tkinter is not installed or not available.")
    sys.exit(1)


def main():
   
    # Step 1: Prompt for input FASTQ
   
    root = tk.Tk()
    root.withdraw()  # Hide main tkinter window

    print("Select the input FASTQ file...")
    input_fastq = filedialog.askopenfilename(
        title="Select input FASTQ",
        filetypes=[("FASTQ files", "*.fq *.fastq"), ("All files", "*.*")]
    )
    if not input_fastq:
        print("No input file selected. Exiting.")
        sys.exit(0)
    print(f"Selected input FASTQ: {input_fastq}")


    # Step 2: Prompt for output FASTQ
 
    print("\nSelect the output FASTQ file (kept reads)...")
    output_fastq = filedialog.asksaveasfilename(
        title="Save trimmed FASTQ as...",
        defaultextension=".fastq",
        filetypes=[("FASTQ files", "*.fq *.fastq"), ("All files", "*.*")]
    )
    if not output_fastq:
        print("No output file selected. Exiting.")
        sys.exit(0)
    print(f"Will save trimmed reads to: {output_fastq}")

    # We'll automatically create a "discarded.fastq" in the same directory:
    discarded_fastq = f"{output_fastq.rsplit('.', 1)[0]}_discarded.fastq"
    print(f"Discarded reads will go to: {discarded_fastq}")


    # Step 3: Ask for Trimming Params

    # Forward primer
    print("\nEnter forward primer (no spaces):")
    fwd_primer_input = input("Forward primer [default: GTCAGCATACTCGAGCTCTG]: ").strip()
    FWD_PRIMER = fwd_primer_input if fwd_primer_input else "GTCAGCATACTCGAGCTCTG"

    print("\nMax spacing before the forward primer (integer)?")
    fwd_spacer_input = input("Default is 10: ").strip()
    MAX_FWD_SPACER = int(fwd_spacer_input) if fwd_spacer_input else 10

    # Reverse primer
    print("\nEnter reverse primer (in forward orientation, as designed):")
    rev_primer_input = input("Reverse primer [default: CTATGCTCGTACTCGCTGTC]: ").strip()
    REV_PRIMER = rev_primer_input if rev_primer_input else "CTATGCTCGTACTCGCTGTC"

    print("\nMax spacing after the reverse primer (integer)?")
    rev_spacer_input = input("Default is 10: ").strip()
    MAX_REV_SPACER = int(rev_spacer_input) if rev_spacer_input else 10

    # Minimum length after trimming
    print("\nEnter the minimum read length after trimming:")
    min_len_input = input("Default is 20: ").strip()
    MIN_LENGTH = int(min_len_input) if min_len_input else 20

    # Maximum length after trimming
    print("\nEnter the maximum read length after trimming:")
    max_len_input = input("Default is 1000: ").strip()
    MAX_LENGTH = int(max_len_input) if max_len_input else 1000


    # Step 4: Build Regex Patterns
    # 1) Forward primer can appear within the first 0..MAX_FWD_SPACER bases
    #    We'll remove everything up to (and including) the primer if matched.
    fwd_pattern = re.compile(rf"^(.{{0,{MAX_FWD_SPACER}}})({FWD_PRIMER})")

    # 2) Reverse primer can appear near the 3' end, within 0..MAX_REV_SPACER of the end
    #    We'll remove everything downstream if matched (including the primer).
    rev_pattern = re.compile(rf"({REV_PRIMER})(.{{0,{MAX_REV_SPACER}}})$")


    # Step 5: Process FASTQ
    num_total = 0
    num_written = 0

    # Counters for reasons of discarding:
    num_discarded_not_trimmed = 0
    num_discarded_too_short = 0
    num_discarded_too_long = 0

    with open(input_fastq, "r") as fin, \
         open(output_fastq, "w") as fout, \
         open(discarded_fastq, "w") as fdisc:

        while True:
            # Read 4 lines (one FASTQ record)
            name_line = fin.readline()
            if not name_line:
                break  # EOF

            seq_line = fin.readline().rstrip("\n")
            plus_line = fin.readline()
            qual_line = fin.readline().rstrip("\n")

            num_total += 1

            original_seq = seq_line
            original_qual = qual_line

          
            # A) Trim forward primer
            match_fwd = fwd_pattern.search(seq_line)
            if match_fwd:
                cut_index_fwd = match_fwd.end(2)  # right after the forward primer
                seq_line = seq_line[cut_index_fwd:]
                qual_line = qual_line[cut_index_fwd:]

            # B) Trim reverse primer
            match_rev = rev_pattern.search(seq_line)
            if match_rev:
                cut_index_rev = match_rev.start(1)  # start of the primer
                seq_line = seq_line[:cut_index_rev]
                qual_line = qual_line[:cut_index_rev]

            # Check if we trimmed anything at all
            was_trimmed = (seq_line != original_seq or qual_line != original_qual)

    
            # C) Filter logic
            #    1) Must be trimmed
            #    2) Must be >= MIN_LENGTH
            #    3) Must be <= MAX_LENGTH
            
            final_len = len(seq_line)

            if not was_trimmed:
                # Discard -> "not trimmed"
                num_discarded_not_trimmed += 1
                # Write read to discarded file
                write_fastq_record(fdisc, name_line, original_seq, plus_line, original_qual)
                continue

            # If was trimmed, check length constraints
            if final_len < MIN_LENGTH:
                num_discarded_too_short += 1
                write_fastq_record(fdisc, name_line, seq_line, plus_line, qual_line)
                continue

            if final_len > MAX_LENGTH:
                num_discarded_too_long += 1
                write_fastq_record(fdisc, name_line, seq_line, plus_line, qual_line)
                continue

            # If passes all filters, write to kept output
            num_written += 1
            write_fastq_record(fout, name_line, seq_line, plus_line, qual_line)


    # Step 6: Print Summary

    num_discarded = num_discarded_not_trimmed + num_discarded_too_short + num_discarded_too_long

    print("\nProcessing complete!")
    print(f"  Total reads processed: {num_total}")
    print(f"  Reads kept (trimmed + within length): {num_written}")
    print(f"  Reads discarded: {num_discarded}")
    print(f"    - Not trimmed: {num_discarded_not_trimmed}")
    print(f"    - Too short (< {MIN_LENGTH}): {num_discarded_too_short}")
    print(f"    - Too long (> {MAX_LENGTH}): {num_discarded_too_long}")
    print(f"\nKept FASTQ: {output_fastq}")
    print(f"Discarded FASTQ: {discarded_fastq}")


def write_fastq_record(handle, name_line, seq, plus_line, qual):
    """Helper function to write 1 FASTQ record to a file handle."""
    handle.write(name_line)
    handle.write(seq + "\n")
    handle.write(plus_line)
    handle.write(qual + "\n")


if __name__ == "__main__":
    main()

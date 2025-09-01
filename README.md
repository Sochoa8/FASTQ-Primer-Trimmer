# FASTQ-Primer-Trimmer
This script provides a way to trim primers from FASTQ files. It was developed as part of my thesis to create a data-cleaning pipeline for Next-Generation Sequencing (NGS) data from SELEX experiments. It serves as a general-purpose adapter trimmer for workflows where primer locations occur after adapters of various nucleotide lengths. 

## Requirements

-   **Python 3.x**
-   **Tkinter**: This library is typically included with standard Python installations on Windows and macOS. On some Linux distributions, it may need to be installed separately.
    ```bash
    # For Debian/Ubuntu
    sudo apt-get install python3-tk
    ```

---

# How to Use 

1. Save the script as cutadapted.py.

2. Open a terminal or command prompt.

3. Navigate to the directory where you saved the file.

4. Run the script using the python command:

```bash
python cutadapted.py
```
The script will then guide you through the interactive steps as described before (file selection via pop-up dialogs, parameter input in the terminal).
---

## Trimming and Filtering Logic 

A read is written to the main output file **only if it meets all three of the following conditions**:

1.  **The read must be trimmed**. At least the forward or reverse primer must be found and removed. Reads where no primers are found are sent to the `discarded.fastq` file.
2.  **The final length must be greater than or equal to `MIN_LENGTH`**. Trimmed reads that are too short are discarded.
3.  **The final length must be less than or equal to `MAX_LENGTH`**. Trimmed reads that are too long are discarded.

This logic ensures that the final output contains only high-confidence sequences that conform to the expected structure of the aptamer library.

### Example Summary Output

After a successful run, you will see a summary like this in your terminal:

```
Processing complete!
  Total reads processed: 100000
  Reads kept (trimmed + within length): 85432
  Reads discarded: 14568
    - Not trimmed: 9876
    - Too short (< 20): 3456
    - Too long (> 1000): 1236

Kept FASTQ: /path/to/your/output_kept.fastq
Discarded FASTQ: /path/to/your/output_discarded.fastq
```

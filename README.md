# Sequence Alignment
This Python script performs DNA sequence alignment using two classic algorithms: Needleman-Wunsch for global alignment and Smith-Waterman for local alignment. The script analyzes predefined DNA sequences under various scoring conditions and visualizes the alignment scoring matrices as heatmaps.

## Key Features & Functionality
The script is built around three core functions that handle the alignment process and output generation.

```needleman_wunsch(seq1, seq2, ...)```
This function implements the Needleman-Wunsch algorithm, a dynamic programming method for global sequence alignment.

- Global Alignment: This approach aligns two sequences across their entire length, attempting to find the best possible match from end to end. It's ideal for comparing sequences that are known to be homologous or highly similar.

- Parameters: It takes two sequences (seq1, seq2) and optional scoring parameters for match_score, mismatch_score, and gap_penalty.

- Output: Returns the two aligned sequences, the full scoring matrix, and the final alignment score.

```smith_waterman(seq1, seq2, ...)```
This function implements the Smith-Waterman algorithm, which is used for local sequence alignment.

- Local Alignment: This method identifies the most similar segments within two sequences. It's highly effective for finding conserved domains or motifs that may be embedded in otherwise unrelated sequences.

- Parameters: Similar to the global alignment function, it takes two sequences and scoring parameters.

- Output: Returns the best aligned subsequences, the full scoring matrix, and the highest local alignment score.

```output(...)```
This utility function processes the results and generates a visual representation of the alignment.

Functionality:

- Prints a summary of the alignment, including the score and a preview of the aligned sequences.

- Generates a heatmap of the scoring matrix using the seaborn and matplotlib libraries, which visually represents the score at each position of the alignment.

- Saves the generated heatmap as a .png file.
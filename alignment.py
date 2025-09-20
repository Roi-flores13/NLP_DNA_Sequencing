import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def needleman_wunsch(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-2):
    n = len(seq1)
    m = len(seq2)
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)
    traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # Initialization
    for i in range(n + 1):
        score_matrix[i, 0] = i * gap_penalty
    for j in range(m + 1):
        score_matrix[0, j] = j * gap_penalty

    # Matrix Fill
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1, j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1, j] + gap_penalty
            insert = score_matrix[i, j - 1] + gap_penalty
            score_matrix[i, j] = max(match, delete, insert)
            
            # Store direction for traceback: 1=diag, 2=up, 3=left
            if score_matrix[i, j] == match:
                traceback_matrix[i, j] = 1
            elif score_matrix[i, j] == delete:
                traceback_matrix[i, j] = 2
            else:
                traceback_matrix[i, j] = 3

    # Traceback to build alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n, m
    while i > 0 or j > 0:
        if traceback_matrix[i, j] == 1: # Diagonal move
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
            
        elif traceback_matrix[i, j] == 2: # Up move (gap in seq2)
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
            
        else: # Left move (gap in seq1)
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
    
    return "".join(reversed(aligned_seq1)), "".join(reversed(aligned_seq2)), score_matrix, score_matrix[n, m]

def smith_waterman(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-2):
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n+1, m+1), dtype=int)
    
    max_score = 0
    max_pos = (0,0)
    
    for i in range(1, n+1):
        for j in range(1, m+1):
            
            match = score_matrix[i-1, j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1, j] + gap_penalty
            insert = score_matrix[i, j-1] + gap_penalty
            
            score_matrix[i, j] = max(0, match, delete, insert)
            
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_pos = (i, j)
                
    aligned_seq1, aligned_seq2 = [], []
    i, j = max_pos
    
    while score_matrix[i, j] > 0:
        current_score = score_matrix[i, j]
        diagonal_score = score_matrix[i-1, j-1]
        up_score = score_matrix[i-1, j]
        left_score = score_matrix[i, j-1]
        
        if current_score == diagonal_score + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
            
        elif current_score == up_score + gap_penalty:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append("-")
            i -= 1
            
        else: 
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j-1])
            j -= 1
            
    return "".join(reversed(aligned_seq1)), "".join(reversed(aligned_seq2)), max_score, score_matrix 

def output(score_nw, aligned_nw1, aligned_nw2, original_seq, variant1_seq, score_matrix_nw, title, plt_nmb):    
    
    print("--- Exercise 1: Needleman-Wunsch Alignment ---")
    print(f"Alignment Score: {score_nw}")
    print(f"Aligned Sequence 1: {aligned_nw1[:100]}...")
    print(f"Aligned Sequence 2: {aligned_nw2[:100]}...")
    print(f"...\n")

    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(score_matrix_nw, annot=False, cmap='viridis', fmt='d')
    plt.title(title)
    plt.xlabel(f"Sequence 2 ({len(variant1_seq)} bp)")
    plt.ylabel(f"Sequence 1 ({len(original_seq)} bp)")
    
    plt.savefig(f"heatmaps/DNA_sequencing_{plt_nmb}.png")
    
    plt.show()
    
plt_number = 1
    
# Sequences
original_seq = """GACTTACGCGCCGTAGCACTTCTGTGATAGCTGCGAGGCGTATTGCTACTTGTACGAGATAGGGTCGACTTT
TCGGAGTCGACAGACACTACGATACT"""
variant1_seq = """AGTATCGTAGTGTCTGTCGACTCCGAAAAGTCGACCCTATCTCGTACAAGTAGCAATACGCCTCGCAGCTAT
CACAGAAGTGCTACGGCGCGTAAGTC"""

original_seq2 = """TTAAATTATATATATACGCGCGCGCGACACACACACACTGATCTATACGCGCGCGCGCGATAGCGATAGCGA
TCGATCGCGCTATATATATA"""
variant1_seq2 = """TATATATAGCGCGATCGATCGCTATCGCTATCGCGCGCGCGTATAGATCAGTGTGTGTGTGTCGCGCGCGCG
TATATATATAATTTAA""" 

original_seq3 = """CTACTGAGTCGTAGCTAGCGAGTCGAGTGCGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT
GCATGCATGCATGCATGC"""
variant1_seq3 = """GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCACGACTCGACTCGCTAGC
TACGACTCAGTAG"""

# Primera Secuencia
#Penalización normal
aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq, variant1_seq)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq, variant1_seq=variant1_seq, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - seq #1 - pen:normal", plt_nmb=plt_number)
plt_number +=1    


aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq, variant1_seq)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq, variant1_seq=variant1_seq, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - se1 #1 - pen:normal", plt_nmb=plt_number)
plt_number +=1

# Penalización match:1, mismatch: -1, gap: -1 

aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq, variant1_seq, match_score=1, gap_penalty=-1)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq, variant1_seq=variant1_seq, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec #1 - pen: 1, -1, -1", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq, variant1_seq, match_score=1, gap_penalty=-1)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq, variant1_seq=variant1_seq, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - sec #1 - pen: 1, -1, -1", plt_nmb=plt_number)
plt_number +=1

# Penalización match:3, mismatch: -2, gap:-1

aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq, variant1_seq, match_score=3, mismatch_score=-2, gap_penalty=-1)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq, variant1_seq=variant1_seq, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec # 1 - pen: 3, -2, -1", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq, variant1_seq, match_score=3, mismatch_score=-2, gap_penalty=-1)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq, variant1_seq=variant1_seq, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - sec # 1 - pen: 3, -2, -1", plt_nmb=plt_number)
plt_number +=1

# Segunda secuencia
# Penalización normal
aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq2, variant1_seq2)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq2, variant1_seq=variant1_seq2, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec #2 - pen:normal", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq2, variant1_seq2)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq2, variant1_seq=variant1_seq2, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment  - sec #2 - pen:normal", plt_nmb=plt_number)
plt_number +=1

# Penalización match:1, mismatch: -1, gap: -1 
aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq2, variant1_seq2, match_score=1, gap_penalty=-1)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq2, variant1_seq=variant1_seq2, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec #2 - pen: 1, -1, -1", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq2, variant1_seq2, match_score=1, gap_penalty=-1)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq2, variant1_seq=variant1_seq2, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - sec #2 - pen: 1, -1, -1", plt_nmb=plt_number)
plt_number +=1

# Penalización match:3, mismatch: -2, gap:-1
aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq2, variant1_seq2, match_score=3, mismatch_score=-2, gap_penalty=-1)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq2, variant1_seq=variant1_seq2, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec #2 - pen: 3, -2, -1", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq2, variant1_seq2, match_score=3, mismatch_score=-2, gap_penalty=-1)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq2, variant1_seq=variant1_seq2, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - sec #2 - pen: 3, -2, -1", plt_nmb=plt_number)
plt_number +=1


#Tercera secuencia
# Penalización normal
aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq3, variant1_seq3)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq3, variant1_seq=variant1_seq3, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec #3 - pen:normal", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq3, variant1_seq3)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq3, variant1_seq=variant1_seq3, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - sec #3 - pen:normal", plt_nmb=plt_number)
plt_number +=1


# Penalización match:1, mismatch: -1, gap: -1 
aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq3, variant1_seq3, match_score=1, gap_penalty=-1)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq3, variant1_seq=variant1_seq3, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec #3 - pen: 1, -1, -1", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq3, variant1_seq3, match_score=1, gap_penalty=-1)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq3, variant1_seq=variant1_seq3, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - sec #3 - pen: 1, -1, -1", plt_nmb=plt_number)
plt_number +=1


# Penalización match:3, mismatch: -2, gap: -1 
aligned_nw1, aligned_nw2, score_matrix_nw, score_nw = needleman_wunsch(original_seq3, variant1_seq3, match_score=3, mismatch_score=-2, gap_penalty=-1)
output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq3, variant1_seq=variant1_seq3, score_matrix_nw=score_matrix_nw, 
       title="Needleman-Wunsch Score Matrix (Global Alignment) - sec #3 - pen: 3, -2, -1", plt_nmb=plt_number)
plt_number +=1

aligned_sw1, aligned_sw2, score_sw, score_matrix_sw = smith_waterman(original_seq3, variant1_seq3, match_score=3, mismatch_score=-2, gap_penalty=-1)

output(score_nw=score_nw, aligned_nw1=aligned_nw1, aligned_nw2=aligned_nw2,
       original_seq=original_seq3, variant1_seq=variant1_seq3, score_matrix_nw=score_matrix_nw,
       title="Smith-Waterman Alignment - sec #3 - pen: 3, -2, -1", plt_nmb=plt_number)
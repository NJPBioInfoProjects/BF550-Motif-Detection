def HammingDistance(s: str, t: str) -> int:
    ham_dist = 0 #Initialize hamming distance value
    for i in range(len(s)): #iterate through s and t at each index and compare
        if s[i] != t[i]: #if s and t do not match at a given index, add 1 to hamming distance count
            ham_dist += 1
    
    return(ham_dist)


CODON_TABLE = {
    'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
    'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C',
    'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
    'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
    'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D',
    'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}

def Translation(s: str) -> str:
    protein = ""  # Start with an empty string
    # Iterate through the RNA string in chunks of 3 (codons)
    for i in range(0, len(s), 3):
        codon = s[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, "")
        if amino_acid == "Stop":
            break  # Stop translation if a stop codon is found
        protein += amino_acid  # Concatenate the amino acid to the protein string
    return protein


def FindingMotif(s: str, t: str) -> list:
    locations = []

    for i in range(len(s) - len(t) + 1): #use +1 to ensure every potential start position of substring is covered
        if s[i:i + len(t)] == t: #see if substring matches
            locations.append(i + 1) #append starting location
    return(locations)

def RNASplicing(dna, introns) -> str:
    dna_str = dna.split('\n')[1] #omit header from fasta
    introns_list = [] #omit header from fasta
    for i in introns:
        introns_list.append(i.split('\n')[1])
    for intron in introns_list:
        dna_str = dna_str.replace(intron, '')  # Remove each intron
    rna = dna_str.replace("T", "U") # Replace T with U to convert to RNA
    
    protein = ""  # Initialize protein sequence
    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]  # Extract each codon
        if codon in CODON_TABLE:
            amino_acid = CODON_TABLE[codon]
            if amino_acid == 'Stop':  # Stop codon encountered, stop translation
                break
            protein += amino_acid  # Add amino acid to the protein string
    return protein

def LongestCommonSubstring(k) -> str:
    #omit headers from fasta
    seqs = []
    for i in k:
        seqs.append(i.split('\n')[1])
    
    #use shortest sequence to prevent indexing error, store length as n
    shortest = min(seqs)
    n = len(shortest)

    for length in range(n, 0, -1):  # From longest to shortest
        for start in range(n - length + 1):  # All possible starting points
            substring = shortest[start:start + length]
            # Check if substring is present in all other sequences
            if all(substring in seq for seq in seqs):
                return substring  # Return as soon as the longest common substring is found
    # If no common substring is found
    return "No common substring found"

def FindingSubsequence(s: str, t: str):
    #omit header from fasta format
    s_split = s.split('\n')[1]
    t_split = t.split('\n')[1]
    indices = []  # To store the 1-indexed positions
    t_index = 0

    for i in range(len(s_split)):
        if t_index >= len(t_split):  # If we have matched all characters in t, break the loop
            break
        if s_split[i] == t_split[t_index]:  # If a character matches
            indices.append(i + 1)  # Store the 1-indexed position
            t_index += 1  # Move to the next character in t
    return(indices)

def LongestCommonSubsequence(s_split: str, t_split: str) -> str:
    #omit header from fasta format
    s_split = s.split('\n')[1]
    t_split = t.split('\n')[1]
    
    m, n = len(s_split), len(t_split)

    # Initialize a 2D lcs_table with all zeroes
    lcs_table = [[0] * (n + 1) for _ in range(m + 1)]

    # Fill the LCS table with the lengths of the LCS
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s_split[i - 1] == t_split[j - 1]:  # If characters match
                lcs_table[i][j] = lcs_table[i - 1][j - 1] + 1
            else:  # Take the maximum of skipping either character
                lcs_table[i][j] = max(lcs_table[i - 1][j], lcs_table[i][j - 1])

    # Reconstruct LCS from the lcs_table
    lcs = []
    i, j = m, n
    while i > 0 and j > 0:
        if s_split[i - 1] == t_split[j - 1]:  # If characters match, add to LCS
            lcs.append(s_split[i - 1])
            i -= 1
            j -= 1
        elif lcs_table[i - 1][j] > lcs_table[i][j - 1]:  # Move up in the table
            i -= 1
        else:  # Move left in the table
            j -= 1

    # Reverse the LCS since it was constructed backwards
    return ''.join(reversed(lcs))


def ShortestCommonSupersequence(s: str, t: str) -> str:
    # Find LCS of s and t using dynamic programming
    m, n = len(s), len(t)
    lcs_table = [[0] * (n + 1) for _ in range(m + 1)]  # Initialize the LCS table

    # Fill the LCS table with the lengths of the LCS
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                lcs_table[i][j] = lcs_table[i - 1][j - 1] + 1
            else:
                lcs_table[i][j] = max(lcs_table[i - 1][j], lcs_table[i][j - 1])

    # Trace back through the LCS table to construct the SCS
    i, j = m, n
    scs = []

    # Trace back to build the Shortest Common Supersequence
    while i > 0 and j > 0:
        if s[i - 1] == t[j - 1]:
            # Add to the SCS and move diagonally if characters match
            scs.append(s[i - 1])
            i -= 1
            j -= 1
        elif lcs_table[i - 1][j] >= lcs_table[i][j - 1]:
            # If coming from above, add s[i-1]
            scs.append(s[i - 1])
            i -= 1
        else:
            # If coming from the left, add t[j-1]
            scs.append(t[j - 1])
            j -= 1

    # Add any remaining characters from s or t
    while i > 0:
        scs.append(s[i - 1])
        i -= 1
    while j > 0:
        scs.append(t[j - 1])
        j -= 1

    # reverse result
    return ''.join(reversed(scs))
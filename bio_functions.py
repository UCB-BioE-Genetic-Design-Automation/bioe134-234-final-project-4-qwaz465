def hairpin_counter(sequence, min_stem=3, min_loop=4, max_loop=9):
    """
    Counts the number of potential hairpin structures in a DNA sequence and returns a simple linear
    representation of the hairpins (stem1(loop)stem2_rc), or None if no hairpins are found.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        min_stem (int): Minimum number of bases in the stem for stable hairpin.
        min_loop (int): Minimum number of bases in the loop.
        max_loop (int): Maximum number of bases in the loop.

    Returns:
        tuple: (int, str or None)
            - The count of potential hairpin structures.
            - A single string showing the detected hairpins in the format 'stem1(loop)stem2_rc', or None if no hairpins are found.
    """
    count = 0
    seq_len = len(sequence)
    hairpin_string = ""

    # Iterate through each base to consider it as a start of a stem
    for i in range(seq_len):
        # Only consider end points that would fit within the loop constraints
        for j in range(i + min_stem + min_loop, min(i + min_stem + max_loop + 1, seq_len)):
            stem1 = sequence[i:i+min_stem]
            stem2 = sequence[j:j+min_stem]

            # Check if the stems are complementary (reverse complement match)
            stem2_rc = reverse_complement(stem2)

            if stem1 == stem2_rc:
                count += 1

                # Extract the loop sequence
                loop = sequence[i+min_stem:j]

                # Create the linear representation (now correctly reversed for output)
                hairpin_representation = f"{stem1}({loop}){stem2}"

                # Append the linear hairpin representation to the string
                hairpin_string += f"Hairpin {count}: {hairpin_representation}\n"

    # Return count and the formatted hairpin string, or None if no hairpins found
    return count, hairpin_string if count > 0 else None

def calculate_edit_distance(s1, s2):
    """
    Compute the edit distance between two strings using a dynamic programming approach based on the Smith-Waterman algorithm for local alignment.

    Parameters:
        s1 (str): The first string to compare.
        s2 (str): The second string to compare.

    Returns:
        int: The edit distance between the two strings, defined as the minimum number of edits (insertions, deletions, or substitutions) required to transform one string into the other.
    """
    s1_len = len(s1)
    s2_len = len(s2)
    dist = [[0] * (s2_len + 1) for _ in range(s1_len + 1)]

    # Initialize distances for transformations involving empty strings
    for i in range(s1_len + 1):
        dist[i][0] = i
    for j in range(s2_len + 1):
        dist[0][j] = j

    # Compute distances
    for i in range(1, s1_len + 1):
        for j in range(1, s2_len + 1):
            if s1[i - 1] == s2[j - 1]:
                dist[i][j] = dist[i - 1][j - 1]
            else:
                dist[i][j] = 1 + min(dist[i - 1][j], dist[i][j - 1], dist[i - 1][j - 1])

    return dist[s1_len][s2_len]

def reverse_complement(sequence):
    """
    Calculates the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The reverse complement of the DNA sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters.
    """
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def translate(sequence):
    """
    Translates a DNA sequence into a protein sequence based on the standard genetic code.

    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The corresponding protein sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters or is not a multiple of three.
    """
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")
    if len(sequence) % 3 != 0:
        raise ValueError("Length of DNA sequence is not a multiple of three, which is required for translation.")

    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        protein += codon_table.get(codon, '_')  # Using '_' for unknown or stop codons
    return protein

if __name__ == "__main__":
    # Example DNA sequence for demonstration
    dna_example = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    try:
        # Perform reverse complement
        rc_result = reverse_complement(dna_example)
        print(f"Reverse Complement of '{dna_example}': {rc_result}")

        # Perform translation
        translation_result = translate(dna_example)
        print(f"Translation of '{dna_example}': {translation_result}")
    except Exception as e:
        print(f"Error: {str(e)}")

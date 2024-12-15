from genome_data_parsing import *
from bio_functions import *
from dataclasses import dataclass
import random
from checkers.forbidden_sequence_checker import ForbiddenSequenceChecker

@dataclass(frozen=True)
class UTROption:
    """
    Encapsulates Ribosome Binding Site (utr) encoding DNAs as modular components (parts) for synthetic biology,
    representing essential sequence elements to facilitate selection algorithms.

    Attributes:
        utr (str): The 5' untranslated region (5' UTR) sequence
        cds (str): The coding sequence of the source gene
        gene_name (str): The name of the source gene
        first_six_aas (str): The precalculated first six amino acids of the source gene's protein sequence
    """
    utr: str
    cds: str
    gene_name: str
    first_six_aas: str

class UTRChooser:
    def __init__(self):
        self.utrOptions = []
        self.kozak_seq = 'aAaAaAATGTCt'
        self.ends = [3,5]
        random.seed(1738)
        self.poly_a_tail_length = 15  # Default length of poly-A tail
    
    def initiate(self):
        self.seq_checker = ForbiddenSequenceChecker()
        self.seq_checker.initiate()
        genbank = 'data/genomic.gbff'
        genes_info = extract_genes_info(genbank)
        file_path = 'data/4932-WHOLE_ORGANISM-integrated.txt'
        top_5_percent_list = proteomics_prune(file_path)

        for locus_tag, abundance in top_5_percent_list:
            if locus_tag in genes_info:
                gene_info = genes_info[locus_tag]

                utr = gene_info['UTR']
                cds = gene_info['CDS']
                gene_name = gene_info['gene']

                # Translate the CDS to get the protein sequence
                try:
                    protein_sequence = translate(str(cds))  # Ensure the CDS is passed as a string
                    # Get the first six amino acids from the protein sequence
                    first_six_aas = protein_sequence[:6]

                    # Create an UTROption instance
                    utr_option = UTROption(
                        utr=str(utr),  # Ensure the UTR is passed as a string
                        cds=str(cds),  # Ensure the CDS is passed as a string
                        gene_name=gene_name,
                        first_six_aas=first_six_aas
                    )

                    # Add the UTROption instance to the class variable list
                    self.utrOptions.append(utr_option)

                except ValueError as e:
                    print(f"Error translating CDS for {locus_tag}: {e}")

    def run(self, cds, end, ignores):
        # in main this should be ran twice to get both 5' and 3' utr
        # end is either 5 or 3 to signify 5' or 3' utr
        # Validate that the CDS is not empty
        if end not in self.ends:
            raise ValueError("Did not specify 3' or 5' end correctly")
        if not cds:
            raise ValueError("CDS sequence cannot be empty.")

        # Validate that the CDS contains only valid nucleotide characters (A, T, C, G)
        if not all(base in 'ATCG' for base in cds):
            raise ValueError("Invalid CDS sequence: contains non-nucleotide characters. Only A, T, C, G are allowed.")

        # Validate that the CDS length is a multiple of 3 for proper translation
        if len(cds) % 3 != 0:
            raise ValueError("CDS sequence length must be a multiple of 3 for valid translation.")

        # Check if there are any utr options to work with
        if not self.utrOptions:
            raise ValueError("No utr options are available to choose from.")

        # Filter out ignored utr options
        valid_utr_options = [utr_option for utr_option in self.utrOptions if utr_option not in ignores]
        if not valid_utr_options:
            raise ValueError("No valid utr options remain after applying the ignore filter.")

        # Initialize the Translate object and translate the first 6 amino acids of the input CDS
        input_first_six_aas = translate(cds[:18])  # 18 bases for 6 amino acids

        best_utr = None
        best_score = float('inf')  # Start with a very high score

        for utr_option in valid_utr_options:
            # check restriction sites

            # if 5' end and doesnt have kozak, skip 
            if end == 5 and not self.ensure_kozak(utr_option, cds):
                continue

            if not self.all_checkers(utr_option):
                continue
            # Calculate the edit distance between input CDS first six AAs and utr option first six AAs
            edit_distance = calculate_edit_distance(input_first_six_aas, utr_option.first_six_aas)

            # Calculate the hairpin count for the UTR + CDS sequence of the utr option
            hairpin_count = hairpin_counter(utr_option.utr + utr_option.cds)[0]

            # Combine the two scores; prioritize the edit distance, then hairpin count
            # doing it this way because if theres a high edit distance then there will be no binding so nothing will happen in the first place
            # however if theres a lot of secondary structure, translation can still occur but much worse, so I am placing much higher priority on
            # edit distance since that more or less says if binding occurs in the first place

            score = (edit_distance * 1000) + hairpin_count  # Giving much higher weight to edit distance

            # If this score is better, update the best utr
            if score < best_score:
                best_score = score
                best_utr = utr_option
        # If selecting a 3' UTR, append the poly-A tail
        if end == 3 and best_utr:
            best_utr.utr += 'A' * self.poly_a_tail_length  # Append poly-A tail to 3' UTR

        return best_utr
    
    def ensure_kozak(self, utr_option, cds):
        #TODO write test for this
        utr = utr_option.utr[-6:]
        seq = utr + cds[:6]
        for x in range(len(seq)):
            # highly conserved
            if self.kozak_seq[x].isupper():
                if self.kozak_seq[x] != seq[x]:
                    return False
            # less conserved, some random chance to pass if not same base
            if self.kozak_seq[x].islower():
                if self.kozak_seq[x].upper() != seq[x]:
                    if random.random() > .3:
                        return False
        return True
    
    def all_checkers(self, utr_option):
        utr = utr_option.utr
        return self.seq_checker.run(utr)

#TODO integrate randomness into kozak checker for lowercase parts
#TODO write benchmarker to see if good?
# Example usage of UTRChooser class

# Create an instance of the UTRChooser
utr_chooser = UTRChooser()

# Initialize the UTRChooser (this loads the genomic data and prepares UTR options)
utr_chooser.initiate()

# Example of using the run method to choose the best UTR for a given CDS
# Assuming we have a CDS sequence, for example, from a gene
cds_sequence = "ATGACCGTACGCTGGAAGGACTGG"  # Sample CDS sequence
end = 5  # We are looking for a 5' UTR
ignores = []  # No UTR options to ignore for this example

# Run the method to choose the best UTR
best_utr_option = utr_chooser.run(cds_sequence, end, ignores)

# Output the chosen UTR Option
if best_utr_option:
    print(f"Best UTR Option for CDS {cds_sequence}:")
    print(f"Gene Name: {best_utr_option.gene_name}")
    print(f"UTR: {best_utr_option.utr}")
    print(f"First Six Amino Acids: {best_utr_option.first_six_aas}")
else:
    print("No suitable UTR option found.")

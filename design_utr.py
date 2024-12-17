from genome_data_parsing import *
from bio_functions import *
from dataclasses import dataclass
import random
from checkers.forbidden_sequence_checker import ForbiddenSequenceChecker

@dataclass(frozen=True)
class UTROption:
    """
    Encapsulates Ribosome Binding Site (UTR) encoding DNAs as modular components for synthetic biology.

    Attributes:
        utr (str): The 5' untranslated region (5' UTR) sequence.
        cds (str): The coding sequence of the source gene.
        gene_name (str): The name of the source gene.
        first_six_aas (str): The precalculated first six amino acids of the source gene's protein sequence.
    """
    utr: str
    cds: str
    gene_name: str
    first_six_aas: str

class UTRChooser:
    """
    A class to select the optimal UTR sequence for a given coding sequence (CDS) based on translation
    compatibility, sequence structure, and scoring criteria.
    
    Attributes:
        kozak_seq (str): A canonical Kozak sequence used for ensuring proper translation initiation.
        ends (list): A list of valid UTR ends (3' or 5').
        poly_a_tail_length (int): The length of the poly-A tail for the 3' UTR.
        utrOptions (list): A list of UTROption instances derived from genomic data.
        seq_checker (ForbiddenSequenceChecker): A sequence checker to validate UTR sequences.
    """
    def __init__(self):
        """
        Initializes the UTRChooser with default settings for Kozak sequence, valid UTR ends, and poly-A tail length.
        """
        self.utrOptions = []
        self.kozak_seq = 'aAaAaAATGTCt'
        self.ends = [3, 5]
        random.seed(1738)
        self.poly_a_tail_length = 15  # Default length of poly-A tail
    
    def initiate(self):
        """
        Loads genomic data and prepares UTR options by parsing gene information and selecting top-performing genes
        based on proteomics data. The UTR options are stored as UTROption instances.
        """
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

                if not utr or not cds:
                    raise ValueError(f"Missing UTR or CDS for gene {gene_name} with locus tag {locus_tag}.")

                if len(cds) % 3 != 0:
                    raise ValueError(f"CDS length for gene {gene_name} is not a multiple of 3.")

                protein_sequence = translate(str(cds))  # Translate CDS to protein sequence
                if len(protein_sequence) < 6:
                    raise ValueError(f"Protein sequence for gene {gene_name} is shorter than 6 amino acids.")

                first_six_aas = protein_sequence[:6]
                
                utr_option = UTROption(
                    utr=str(utr),
                    cds=str(cds),
                    gene_name=gene_name,
                    first_six_aas=first_six_aas
                )
                self.utrOptions.append(utr_option)

    def run(self, cds, end, ignores):
        """
        Selects the best UTR option for a given coding sequence (CDS) based on scoring criteria.

        Parameters:
            cds (str): The coding sequence for which the UTR is being chosen.
            end (int): Specifies 5' or 3' UTR (use 5 or 3).
            ignores (set): A set of UTROption instances to exclude from selection.

        Returns:
            UTROption: The best UTR option for the given CDS.
        """
        if end not in self.ends:
            raise ValueError("End must be 3 or 5 to signify 3' or 5' UTR.")
        if not cds:
            raise ValueError("CDS sequence cannot be empty.")
        if not all(base in 'ATCG' for base in cds):
            raise ValueError("CDS sequence contains invalid characters. Only A, T, C, G are allowed.")
        if len(cds) % 3 != 0:
            raise ValueError("CDS sequence length must be a multiple of 3.")
        if len(cds) < 18:
            raise ValueError("CDS sequence is too short to translate the first six amino acids.")
        if not self.utrOptions:
            raise ValueError("No UTR options are available to choose from.")
        
        valid_utr_options = [utr_option for utr_option in self.utrOptions if utr_option not in ignores]
        if not valid_utr_options:
            raise ValueError("No valid UTR options remain after applying the ignore filter.")

        input_first_six_aas = translate(cds[:18])  # First 6 amino acids from the CDS
        best_utr = None
        best_score = float('inf')

        kozak_compliant_found = False

        for utr_option in valid_utr_options:
            if end == 5:
                is_kozak_compliant = self.ensure_kozak(utr_option, cds)
                kozak_compliant_found = kozak_compliant_found or is_kozak_compliant
                if not is_kozak_compliant:
                    continue  # Skip options that do not meet Kozak sequence criteria
            if not self.all_checkers(utr_option):
                continue  # Skip options failing forbidden sequence checks

            edit_distance = calculate_edit_distance(input_first_six_aas, utr_option.first_six_aas)
            hairpin_count = hairpin_counter(utr_option.utr + utr_option.cds)[0]

            # Weighted scoring: edit distance has higher priority
            score = (edit_distance * 1000) + hairpin_count

            if score < best_score:
                best_score = score
                best_utr = utr_option

        if end == 5 and not kozak_compliant_found:
            raise ValueError("No Kozak-compliant UTR options found.")

        if end == 3 and best_utr:
        # Append poly-A tail for 3' UTR by creating a new UTROption
            best_utr = UTROption(
            utr=best_utr.utr + 'A' * self.poly_a_tail_length,  # Append poly-A tail
            cds=best_utr.cds,
            gene_name=best_utr.gene_name,
            first_six_aas=best_utr.first_six_aas
        )

        return best_utr
    
    def ensure_kozak(self, utr_option, cds):
        """
        Ensures the UTR sequence ends with a Kozak-like sequence before the start codon.

        Parameters:
            utr_option (UTROption): The UTR option to validate.
            cds (str): The coding sequence.

        Returns:
            bool: True if the UTR contains a valid Kozak sequence, False otherwise.
        """
        utr = utr_option.utr[-6:]
        seq = utr + cds[:6]

        for x in range(len(seq)):
            if self.kozak_seq[x].isupper() and self.kozak_seq[x] != seq[x]:
                return False  # Highly conserved region mismatch
            if self.kozak_seq[x].islower() and self.kozak_seq[x].upper() != seq[x]:
                if random.random() > 0.3:  # Less conserved region with chance to pass
                    return False
        return True
    
    def all_checkers(self, utr_option):
        """
        Validates the UTR sequence against forbidden sequence checks.

        Parameters:
            utr_option (UTROption): The UTR option to validate.

        Returns:
            bool: True if the sequence passes all checks, False otherwise.
        """
        return self.seq_checker.run(utr_option.utr)

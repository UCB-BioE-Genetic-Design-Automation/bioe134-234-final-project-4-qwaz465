from Bio.SeqUtils import MeltingTemp as mt
from bio_functions import reverse_complement

class PrimerDesigner:
    """
    A class to design primers for molecular cloning experiments using 
    either Golden Gate cloning or Gibson Assembly.
    """
    
    def __init__(self):
        """
        Initialize the PrimerDesigner class with default settings.
        Supported cloning methods are 'Golden Gate' and 'Gibson'.
        """
        self.enzyme_dict = {}
        self.overhangs = tuple()
        self.methods = []

    def initiate(self):
        """
        Initialize the default parameters for the primer design process:
        - Target melting temperature (Tm)
        - Primer length
        - Restriction enzyme recognition sites (Golden Gate)
        - Homology region length (Gibson)
        """
        self.tm_target = 60
        self.primer_length = 20
        self.enzyme_dict = {
            "BsaI": "GGTCTC",
            "BsmBI": "CGTCTC",
            "BbsI": "GAAGAC",
            "Esp3I": "GCGGCC",
            "AarI": "CACCTGC",
            "SapI": "GCTCTTC"
        }
        self.overhangs = ("AGCT", "TCGA")  # Default overhang sequences
        self.homology_length = 30           # Default homology length for Gibson Assembly
        self.methods = ['Gibson', 'Golden Gate']

    def run(self, cds, utr5, utr3, enzyme=None, method="Golden Gate"):
        """
        Generate primers for a given sequence using the selected cloning method.

        Parameters:
        - cds: str, Coding sequence.
        - utr5: str, 5' UTR sequence.
        - utr3: str, 3' UTR sequence.
        - enzyme: str, Restriction enzyme name (only for Golden Gate).
        - method: str, Cloning method ('Golden Gate' or 'Gibson').

        Returns:
        - dict: Forward and reverse primers.
        """
        if method not in self.methods:
            raise ValueError(f"Unsupported method, choose from {self.methods}")
        
        full_sequence = utr5 + cds + utr3

        if method == "Golden Gate":
            if enzyme not in self.enzyme_dict:
                raise ValueError(f"Unsupported enzyme: {enzyme}. Choose from: {list(self.enzyme_dict.keys())}")
            cloning_site = self.enzyme_dict[enzyme]

            # Forward primer: Overhang + Cloning site + Start of the sequence
            forward_primer_core = self._adjust_primer_length(full_sequence[:self.primer_length])
            forward_primer = self.overhangs[0] + cloning_site + forward_primer_core

            # Reverse primer: Overhang + Cloning site + Reverse complement of the end of the sequence
            reverse_primer_core = self._adjust_primer_length(full_sequence[-self.primer_length:])
            reverse_primer_core_rc = reverse_complement(reverse_primer_core)
            reverse_primer = self.overhangs[1] + cloning_site + reverse_primer_core_rc

        elif method == "Gibson":
            # Forward primer: Homology region + Start of the sequence
            forward_primer_core = self._adjust_primer_length(full_sequence[:self.primer_length])
            forward_primer = self._add_homology_region("upstream", forward_primer_core)

            # Reverse primer: Homology region + Reverse complement of the end of the sequence
            reverse_primer_core = self._adjust_primer_length(full_sequence[-self.primer_length:])
            reverse_primer_core_rc = reverse_complement(reverse_primer_core)
            reverse_primer = self._add_homology_region("downstream", reverse_primer_core_rc)
    
        return {
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer
        }

    def _adjust_primer_length(self, sequence):
        """
        Adjust the primer length to meet the target melting temperature (Tm).

        Parameters:
        - sequence: str, Sequence to trim and evaluate.

        Returns:
        - str: Sequence adjusted to meet the Tm target.
        """
        length = self.primer_length
        while length < len(sequence):
            tm = mt.Tm_NN(sequence[:length])
            if tm >= self.tm_target:
                return sequence[:length]
            length += 1
        return sequence

    def _add_homology_region(self, region_type, primer_core):
        """
        Add homology regions to primers for Gibson Assembly.

        Parameters:
        - region_type: str, Either 'upstream' or 'downstream' to indicate location.
        - primer_core: str, The core primer sequence.

        Returns:
        - str: Primer with added homology region.
        """
        if region_type == "upstream":
            homology = "A" * self.homology_length  # Replace with real upstream sequence
        elif region_type == "downstream":
            homology = "T" * self.homology_length  # Replace with real downstream sequence
        else:
            raise ValueError("Invalid region type. Use 'upstream' or 'downstream'.")
        return homology + primer_core

# Example usage:
# Golden Gate
golden_gate_designer = PrimerDesigner()
golden_gate_designer.initiate()
primers_golden_gate = golden_gate_designer.run(
    cds="ATGACCTGACTGA",
    utr5="TTTAAA",
    utr3="TTTCCC",
    enzyme="BsaI",
    method="Golden Gate"
)
print("Golden Gate Primers:", primers_golden_gate)

# Gibson Assembly
gibson_designer = PrimerDesigner()
gibson_designer.initiate()
primers_gibson = gibson_designer.run(
    cds="ATGACCTGACTGA",
    utr5="TTTAAA",
    utr3="TTTCCC",
    method="Gibson"
)
print("Gibson Assembly Primers:", primers_gibson)

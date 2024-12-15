from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from bio_functions import *
class PrimerDesigner:
    # cut out cds's from utr genes and only take utrs from them
    # also use specifically 1 cloning tech for nowwwww
    #currently for golden gate
    def __init__(self):
        self.enzyme_dict = {}
        self.overhangs = tuple()

    def initiate(self):
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
        self.overhangs = ("AGCT", "TCGA")
    def run(self, cds, utr5, utr3, enzyme):
        if enzyme not in self.enzyme_dict:
            raise ValueError(f"Unsupported enzyme: {enzyme}. Please choose from: {list(self.enzyme_dict.keys())}")
        # Generate the forward primer
        # Prepare full sequence
        full_sequence = utr5 + cds + utr3
        cloning_site = self.enzyme_dict[enzyme]

        # Forward primer: Overhang + Cloning site + Start of the sequence
        forward_primer_core = self._adjust_primer_length(full_sequence[:self.primer_length])
        forward_primer = self.overhangs[0] + cloning_site + forward_primer_core

        # Reverse primer: Overhang + Cloning site + Reverse complement of the end of the sequence
        reverse_primer_core = self._adjust_primer_length(full_sequence[-self.primer_length:])
        reverse_primer_core_rc = reverse_complement(reverse_primer_core)
        reverse_primer = self.overhangs[1] + cloning_site + reverse_primer_core_rc

        return {
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer
        }
    
    def _adjust_primer_length(self, sequence):
        """Adjust the primer length to match the target Tm."""
        length = self.primer_length
        while length < len(sequence):
            tm = mt.Tm_NN(sequence[:length])
            if tm >= self.tm_target:
                return sequence[:length]
            length += 1
        return sequence
    
    # Example usage:
designer = PrimerDesigner(
    cds="ATGACCTGACTGA",
    utr5="TTTAAA",
    utr3="TTTCCC",
    tm_target=60,
    primer_length=20,
    enzyme="BsaI",
    overhangs=("AGCT", "TCGA")
)
primers = designer.run()
print(primers)
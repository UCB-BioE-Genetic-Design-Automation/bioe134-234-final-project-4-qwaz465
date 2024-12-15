from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
class PrimerDesigner:
    # cut out cds's from utr genes and only take utrs from them
    # also uae specifically 1 cloning tech for nowwwww
    def initiate(self):
        self.tm_target = 60
        self.primer_length = 20
    def run(self, cds, utr5, utr3, cloning_site):
        # Generate the forward primer
        full_sequence = utr5 + cds + utr3
        forward_primer_core = full_sequence[:self.primer_length]
        while mt.Tm_NN(forward_primer_core) < self.tm_target and len(forward_primer_core) < len(full_sequence):
            forward_primer_core += full_sequence[len(forward_primer_core)]
        forward_primer = forward_primer_core

        # Add cloning site 
        forward_primer = cloning_site + forward_primer

        # Generate the reverse primer (complement of the reverse strand)
        reverse_primer_core = Seq(full_sequence[-self.primer_length:]).reverse_complement()
        while mt.Tm_NN(reverse_primer_core) < self.tm_target and len(reverse_primer_core) < len(full_sequence):
            reverse_primer_core = Seq(full_sequence[-(len(reverse_primer_core) + 1):]).reverse_complement()
        reverse_primer = str(reverse_primer_core)

        # Add cloning site
        reverse_primer = cloning_site + reverse_primer

        return forward_primer, reverse_primer
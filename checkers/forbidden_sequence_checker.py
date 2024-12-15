from bio_functions import *
class ForbiddenSequenceChecker:
    def __init__(self):
        self.forbidden = []

    def initiate(self):
        # Populate forbidden sequences
        self.forbidden = [
    "AAAAAAAA",  # poly(A)
    "TTTTTTTT",  # poly(T)
    "CCCCCCCC",  # poly(C)
    "GGGGGGGG",  # poly(G)
    "ATATATAT",  # poly(AT)
    "GUAAGU",    # Splice donor site (cryptic)
    "CAG",       # Splice acceptor site (cryptic)
    "GGGGGG",    # G-quadruplex motif
    "ATATATATATAT",  # AT-rich region (over 75% AT)
    "TTTTTTTTTT",  # Repetitive T's
    "AATAAA",    # Polyadenylation signal (AATAAA)
    "TATAAA",    # Polyadenylation signal (TATAAA)
    "TATAAA",    # TATA box (promoter-like)
    "AATAAA",    # Terminator-like sequence
    "GGAGGGGAGAG", # Ty1 Transposon sequence
    "TGAGGGGG",  # LTR sequence
    "AGGAGG",    # Shine-Dalgarno-like sequence
    "TAG",       # Stop codon
    "TAA",       # Stop codon
    "TGA",       # Stop codon
    "CAATTG",    # MfeI
    "GAATTC",    # EcoRI
    "GGATCC",    # BamHI
    "AGATCT",    # BglII
    "ACTAGT",    # SpeI
    "TCTAGA",    # XbaI
    "GGTCTC",    # BsaI
    "CGTCTC",    # BsmBI
    "CACCTGC",   # AarI
    "CTGCAG",    # PstI
    "CTCGAG",    # XhoI
    "GCGGCCGC",  # NotI
    "AAGCTT",    # HindIII
]


    def run(self, dnaseq):
        # Use the reverse_complement function from seq_utils
        rc = reverse_complement(dnaseq)
        combined = (dnaseq + "x" + rc).upper()

        for site in self.forbidden:
            if site in combined:
                return False

        return True
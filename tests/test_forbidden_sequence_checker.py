import pytest
from checkers.forbidden_sequence_checker import ForbiddenSequenceChecker

@pytest.fixture
def checker():
    checker = ForbiddenSequenceChecker()
    checker.initiate()
    return checker

def test_forbidden_sites(checker):
    # Sequences with forbidden sites
    forbidden_seqs = [
        "TTGACAATTgaattcCGAACTAGTATAAT",
        "TTTTTTTTTTGTCGAGAAATTTATAAT",
        "TTGACATTACCGTCTCGAGCGCCTATAAT",
        "TTGACAGCGAACGCTTCAGACTGCAGAT",
        "TTGACTGCAGTTGTAACTTATATAAT",
        "TTGACAAGAAGCGGCCGCTCAATTATAAT",
        "TTGACATTATGACACCTGCTTATTATAAT",
        "TTGACACTCGAGCACAGGCTCTATAAT",
        "TTGACATCGGGGGGGGTTTTACCATGGTCGTTATAAT",
        "TTGACAAAGTCGATTTTTTTTTTCCTTCGATTATAAT",
        "TTGACACGGTCTCATTCACTAGGTTATAAT",
        "TTGACAGTCTAGAGTCTGAACAAGGAGATCTTAAT",
    ]

    print(">> Testing problematic sequences (expected False)")
    for seq in forbidden_seqs:
        result = checker.run(seq.upper())
        print(f"result: {result} on {seq}")
        assert result == False

def test_allowed_sites(checker):
    # Sequences without forbidden sites
    allowed_seqs = [
        "AAACTGTAATCCACCACAAGTCAAGCCAT",
        "GCCTCTCTGAGGACGCCGTATGAATTAATA",
        "GTAAACTTTGCGCGGGTTCACTGCGATCC",
        "TTCAGTCTCGTCCAAGGGCACAATCGAAT",
        "ATCCCCCGAAGTTTAGCAGGTCGTGAGGT",
        "TCATGGAGGCTCTCGTTCATCCCGTGGGA",
        "ATCAAGGCTTCGCCTTGATAAAGCACCCCG",
        "TCGGGTGTAGCAGAGAAGACGCCTACTGA",
        "TTGTGCGATCCCTCCACCTCAGCTAAGGT",
        "GCTACCAATATTTAGTTTTTTAGCCTTGC",
        "ACAGACATCCTACTTAGATTGCCACGCAT",
        "GTTCCGCTGGCGATCCATCGTTGGCGGCCG",
    ]

    print("\n>> Testing random sequences (expected True)")
    for seq in allowed_seqs:
        result = checker.run(seq)
        print(f"result: {result} on {seq}")
        assert result == True

def test_repeated_forbidden_sequence(checker):
    # Test case with repeated instances of a forbidden sequence (e.g., "AATAAA")
    seq_with_repeated_forbidden = "AATAAA" * 5  # "AATAAA" repeated
    result = checker.run(seq_with_repeated_forbidden)
    print(f"result: {result} on {seq_with_repeated_forbidden}")
    assert result == False

def test_non_forbidden_sequence_with_similar_substrings(checker):
    # Test a sequence similar to a forbidden sequence but without an exact match
    similar_seq = "GCGGGCGCTCGTTCGAGTATAAT"
    result = checker.run(similar_seq)
    print(f"result: {result} on {similar_seq}")
    assert result == True

def test_check_forbidden_with_special_characters(checker):
    # Check a sequence with special characters; they should not interfere with the check
    special_char_seq = "AAATAA$%&@#AATAA"
    with pytest.raises(ValueError):
        checker.run(special_char_seq)

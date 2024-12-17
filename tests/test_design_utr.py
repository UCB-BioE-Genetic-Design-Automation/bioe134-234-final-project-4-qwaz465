import pytest
from design_utr import UTRChooser, UTROption

@pytest.fixture
def utr_chooser():
    chooser = UTRChooser()
    chooser.initiate()
    return chooser

def test_invalid_end(utr_chooser):
    # Test invalid end value (should be 3 or 5)
    with pytest.raises(ValueError, match="End must be 3 or 5 to signify 3' or 5' UTR."):
        utr_chooser.run("ATGCGTA", 7, set())

def test_empty_cds(utr_chooser):
    # Test empty CDS sequence
    with pytest.raises(ValueError, match="CDS sequence cannot be empty."):
        utr_chooser.run("", 5, set())

def test_invalid_cds_characters(utr_chooser):
    # Test CDS sequence with invalid characters
    with pytest.raises(ValueError, match="CDS sequence contains invalid characters. Only A, T, C, G are allowed."):
        utr_chooser.run("ATGXXGTA", 5, set())

def test_cds_not_multiple_of_three(utr_chooser):
    # Test CDS length that is not a multiple of 3
    with pytest.raises(ValueError, match="CDS sequence length must be a multiple of 3."):
        utr_chooser.run("ATGAT", 5, set())

def test_cds_too_short(utr_chooser):
    # Test CDS that is too short
    with pytest.raises(ValueError, match="CDS sequence is too short to translate the first six amino acids."):
        utr_chooser.run("ATG", 5, set())

def test_no_utr_options(utr_chooser):
    # Test when there are no UTR options available
    utr_chooser.utrOptions = []  # Empty the UTR options list
    with pytest.raises(ValueError, match="No UTR options are available to choose from."):
        utr_chooser.run("ATGTCTGCGGGCGCTCGTTCGAGTATAATC", 5, set())

def test_no_kozak_compliant_utr(utr_chooser):
    # Test when no UTR options are Kozak compliant
    invalid_utr = UTROption(utr="GGGGAAGG", cds="ATGCATG", gene_name="GeneX", first_six_aas="MALQ")
    utr_chooser.utrOptions = [invalid_utr]
    with pytest.raises(ValueError, match="No Kozak-compliant UTR options found."):
        utr_chooser.run("ATGTGTGCGGGCGCTCGTTCGAGTATAATC", 5, set())

def test_valid_utr_output(utr_chooser):
    # Test valid UTR selection
    valid_utr = UTROption(utr="ACGGACGGTCCACCTAAAAAA", cds="ATGCATG", gene_name="GeneX", first_six_aas="MALQ")
    utr_chooser.utrOptions = [valid_utr]
    
    result = utr_chooser.run("ATGTCTGCGGGCGCTCGTTCGAGTATAATC", 5, set())
    
    # Ensure that the correct output is returned
    assert result == valid_utr
    assert result.utr == "ACGGACGGTCCACCTAAAAAA"  # Ensure the UTR is correctly assigned
    
def test_poly_a_tail_for_3_end(utr_chooser):
    # Test that poly-A tail is added to the 3' UTR
    valid_utr = UTROption(utr="ACGGACGGTCCACCTAAAAAA", cds="ATGCATG", gene_name="GeneX", first_six_aas="MALQ")
    utr_chooser.utrOptions = [valid_utr]
    
    result = utr_chooser.run("ATGTCTGCGGGCGCTCGTTCGAGTATAATC", 3, set())
    
    # Ensure the poly-A tail is appended
    assert result.utr == "ACGGACGGTCCACCTAAAAAA" + "A" * utr_chooser.poly_a_tail_length  # Check poly-A tail is added correctly

def test_all_checkers_pass(utr_chooser):
    # Test all checkers pass (no forbidden sequences)
    valid_utr = UTROption(utr="ACGGACGGTCCACCTAAAAAA", cds="ATGCATG", gene_name="GeneX", first_six_aas="MALQ")
    utr_chooser.utrOptions = [valid_utr]
    
    result = utr_chooser.run("ATGTCTGCGGGCGCTCGTTCGAGTATAATC", 5, set())
    
    # Ensure the UTR passes all checks and is returned as valid
    assert result == valid_utr

def test_forbidden_sequence_fail(utr_chooser):
    # Test forbidden sequence checker failure
    forbidden_utr = UTROption(utr="AAAAAAAA", cds="ATGCATG", gene_name="GeneX", first_six_aas="MALQ")
    utr_chooser.utrOptions = [forbidden_utr] + utr_chooser.utrOptions
    
    result = utr_chooser.run("ATGTCTGCGGGCGCTCGTTCGAGTATAATC", 5, set())
    
    # Assert the UTR is excluded due to forbidden sequence
    assert result != forbidden_utr

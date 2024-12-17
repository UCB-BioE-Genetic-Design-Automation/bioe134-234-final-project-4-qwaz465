import pytest
from design_utr import *

@pytest.fixture
def utr_chooser():
    chooser = UTRChooser()
    chooser.initiate()
    return chooser

def test_kozak_no(utr_chooser):
    # Test where the UTR and CDS don't match the Kozak sequence
    utr_option_fail = UTROption(
        utr="AAACCTA",  # UTR doesn't match, last 'A' is incorrect for Kozak 'A'
        cds="ATGGCACG",  # CDS doesn't match for the required sequence (first 6 AAs won't match)
        gene_name="GeneFail",
        first_six_aas="MAAR"
    )
    cds = "ATGGCACG"  # Ensure the CDS sequence doesn't match the Kozak sequence
    assert utr_chooser.ensure_kozak(utr_option_fail, cds) == False

def test_kozak_yes(utr_chooser):
    # Test where the UTR and CDS match the Kozak sequence
    utr_option_pass = UTROption(
        utr="AAAAAAA",  # UTR matches Kozak sequence: 'A' at the uppercase positions, 'a' at the flexible ones
        cds="ATGTCCG",  # CDS matches Kozak sequence: 'T' and 'C' at positions 4 and 5
        gene_name="GenePass",
        first_six_aas="MAAR"
    )
    cds = "ATGTCTG"  # Ensure the CDS sequence matches the required pattern in the Kozak sequence
    assert utr_chooser.ensure_kozak(utr_option_pass, cds) == True

def test_kozak_with_random_chance(utr_chooser):
    # Test with a UTR that has some mismatches but passes due to random chance
    utr_option_rand = UTROption(
        utr="AAAAAAA",  # UTR is the same as the Kozak sequence
        cds="ATGGCACG",  # CDS does not exactly match the Kozak sequence
        gene_name="GeneRand",
        first_six_aas="MAAR"
    )
    cds = "ATGGCACG"  # Ensure the CDS doesn't fully match but might pass due to random chance
    result = utr_chooser.ensure_kozak(utr_option_rand, cds)
    # Check if result is either True or False based on random chance
    # with seed this always fails
    assert result == False
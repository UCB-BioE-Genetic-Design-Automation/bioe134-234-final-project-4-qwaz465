import pytest
from design_utr import *

@pytest.fixture
def utr_chooser():
    chooser = UTRChooser()
    chooser.initiate()
    return chooser

def test_kozak_no(utr_chooser):
    utr_option_pass = UTROption(
    utr="AAACCCA",
    cds="ATGGCCGCT",
    gene_name="GeneFail",
    first_six_aas="MAAR"
)
    cds = "ATGGCCGCT"
    assert utr_chooser.ensure_kozak(utr_option_pass, cds) == False

def test_kozak_yes(utr_chooser):
    utr_option_pass = UTROption(
    utr="AAACACA",
    cds="ATGTCCGCT",
    gene_name="GenePass",
    first_six_aas="MAAR"
)
    cds = "ATGTCCGCT"
    assert utr_chooser.ensure_kozak(utr_option_pass, cds) == True
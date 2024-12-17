import pytest
from design_primer import PrimerDesigner  # Assuming PrimerDesigner class is in primer_designer.py

@pytest.fixture
def primer_designer():
    # Fixture to initialize the PrimerDesigner
    designer = PrimerDesigner()
    designer.initiate()
    return designer

def test_invalid_method(primer_designer):
    # Test that ValueError is raised if an unsupported method is used
    with pytest.raises(ValueError, match="Unsupported method, choose from"):
        primer_designer.run(
            cds="ATGACCTGACTGA",
            utr5="TTTAAA",
            utr3="TTTCCC",
            method="UnsupportedMethod"
        )

def test_invalid_enzyme(primer_designer):
    # Test that ValueError is raised if an invalid enzyme is used for Golden Gate method
    with pytest.raises(ValueError, match="Unsupported enzyme:"):
        primer_designer.run(
            cds="ATGACCTGACTGA",
            utr5="TTTAAA",
            utr3="TTTCCC",
            enzyme="InvalidEnzyme",
            method="Golden Gate"
        )

def test_valid_golden_gate_primers(primer_designer):
    # Test that valid primers are generated for Golden Gate method
    primers = primer_designer.run(
        cds="ATGACCTGACTGA",
        utr5="TTTAAA",
        utr3="TTTCCC",
        enzyme="BsaI",
        method="Golden Gate"
    )

    # Check that both forward and reverse primers are returned
    assert "forward_primer" in primers
    assert "reverse_primer" in primers

    # Check that primers have the expected structure (e.g., overhang + cloning site + core)
    assert primers["forward_primer"].startswith("AGCTGGTCTC")
    assert primers["reverse_primer"].startswith("TCGAGGTCTC")

def test_valid_gibson_primers(primer_designer):
    # Test that valid primers are generated for Gibson method
    primers = primer_designer.run(
        cds="ATGACCTGACTGA",
        utr5="TTTAAA",
        utr3="TTTCCC",
        method="Gibson"
    )

    # Check that both forward and reverse primers are returned
    assert "forward_primer" in primers
    assert "reverse_primer" in primers

    # Check that primers include homology regions
    assert primers["forward_primer"].startswith("A" * 30)  # Homology region for upstream
    assert primers["reverse_primer"].startswith("T" * 30)  # Homology region for downstream

def test_invalid_homology_region(primer_designer):
    # Test that ValueError is raised if an invalid region type is used for Gibson
    with pytest.raises(ValueError, match="Invalid region type. Use 'upstream' or 'downstream'."):
        primer_designer._add_homology_region(
            region_type="invalid_region",
            primer_core="ATGACCTGACTGA"
        )

def test_primer_length_adjustment(primer_designer):
    # Test that the primers are adjusted to meet the target Tm
    primers = primer_designer.run(
        cds="ATGACCTGACTGA",
        utr5="TTTAAA",
        utr3="TTTCCC",
        enzyme="BsaI",
        method="Golden Gate"
    )

    # Check that the length of primers is adjusted correctly
    forward_primer = primers["forward_primer"]
    reverse_primer = primers["reverse_primer"]

    # Ensure the primer lengths are as expected (20 bases by default)
    assert len(forward_primer) >= 20
    assert len(reverse_primer) >= 20
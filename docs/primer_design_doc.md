# Yeast Primer Design: `design_primer.py`

The `design_primer.py` script is part of a Python-based pipeline for designing plasmid prototypes for yeast. It generates primers compatible with **Golden Gate Assembly** and **Gibson Assembly**, tailoring the design logic to the requirements of each method.

## Overview

### Key Features:
- **Golden Gate Assembly:**  
  - Restriction enzyme sites are added from a predefined dictionary (`enzyme_dict`).
  - Specific overhanging nucleotides are appended to the primer sequence.
  
- **Gibson Assembly:**  
  - Homology regions are incorporated using `_add_homology_region()` instead of restriction sites.

- **Melting Temperature Adjustment:**  
  - Primer lengths are dynamically adjusted using `_adjust_primer_length()` to meet the target melting temperature (Tm).

- **Sequence Inputs:**  
  - The input sequence includes two untranslated regions (UTRs) and a coding sequence (CDS).

---

## Primer Design Logic

### 1. **Golden Gate Assembly**
- Restriction enzyme recognition sites are added to the primers.
- Overhang sequences (`self.overhangs`) are prepended for compatibility.
- Primers are stitched together in the following format:  
  **Overhang + Restriction Site + Adjusted Primer Core.**

### 2. **Gibson Assembly**
- Homology regions are added upstream and downstream of the primer core using `_add_homology_region()`.
- Homology region lengths are set via `self.homology_length`.

---

## Key Functions

### `_adjust_primer_length()`
- Dynamically adjusts primer length to achieve the desired melting temperature (Tm).
- Uses Biopython's `MeltingTemp` module (`Tm_NN`) for Tm calculations.

### `_add_homology_region()`
- Adds homology regions to primers for Gibson Assembly.  
- Prepares homology sequences depending on the region type:
  - `"upstream"`: Prepends a default or custom upstream homology sequence.
  - `"downstream"`: Prepends a default or custom downstream homology sequence.

---

## Workflow

1. **Input:**  
   Provide the following components:  
   - **CDS (coding sequence)**  
   - **UTR5 (5' untranslated region)**  
   - **UTR3 (3' untranslated region)**  
   - **Method** (`"Golden Gate"` or `"Gibson"`)  
   - **Enzyme** (optional, required for Golden Gate).

2. **Processing:**  
   - Combine UTRs and CDS into a full sequence.  
   - Generate primers based on the specified method.

3. **Output:**  
   A dictionary containing:  
   - `forward_primer`: The forward primer sequence.  
   - `reverse_primer`: The reverse primer sequence.  

---

---

#### **Error Handling**

The `PrimerDesigner` class includes error handling mechanisms to ensure proper use of the tool and prevent common issues during primer design. Below are some key error scenarios and how they are handled:

1. **Unsupported Cloning Method**:
   - **Error**: `ValueError: Unsupported method, choose from ['Gibson', 'Golden Gate']`
   - **Cause**: The cloning method passed to the `run` method is not supported. The valid methods are "Golden Gate" and "Gibson".
   - **Resolution**: Ensure that the method provided is either "Golden Gate" or "Gibson". If you're unsure about the method, check the available options in the `methods` attribute.

2. **Invalid Enzyme for Golden Gate Cloning**:
   - **Error**: `ValueError: Unsupported enzyme: BsaI. Choose from: ['BsaI', 'BsmBI', 'BbsI', 'Esp3I', 'AarI', 'SapI']`
   - **Cause**: The enzyme passed for Golden Gate cloning is not available in the `enzyme_dict`.
   - **Resolution**: When using Golden Gate cloning, ensure the enzyme is one of the recognized restriction enzymes (e.g., "BsaI", "BsmBI", etc.). Check the `enzyme_dict` for the list of supported enzymes.

3. **Invalid Cloning Method for Gibson**:
   - **Error**: `ValueError: Invalid region type. Use 'upstream' or 'downstream'.`
   - **Cause**: The `region_type` parameter passed to the `_add_homology_region` method is neither "upstream" nor "downstream".
   - **Resolution**: When using Gibson assembly, ensure that the `region_type` is set correctly to either "upstream" or "downstream". This specifies whether the homology region is upstream or downstream of the primer sequence.

By handling these errors gracefully, the `PrimerDesigner` class provides clear feedback on what went wrong during the primer design process, allowing users to troubleshoot and fix issues effectively.

## Example Usage

### Golden Gate Assembly:
```python
designer = PrimerDesigner()
designer.initiate()
primers = designer.run(
    cds="ATGACCTGACTGA",
    utr5="TTTAAA",
    utr3="TTTCCC",
    enzyme="BsaI",
    method="Golden Gate"
)
print("Golden Gate Primers:", primers)
```

**Output:**
```python
Golden Gate Primers: {
    'forward_primer': 'AGCTGGTCTCTTTAAATGACCTGACTGA',
    'reverse_primer': 'TCGAGGTCTCGGGAAAGTCAGGTCAT'
}
```

---

### Gibson Assembly:
```python
designer = PrimerDesigner()
designer.initiate()
primers = designer.run(
    cds="ATGACCTGACTGA",
    utr5="TTTAAA",
    utr3="TTTCCC",
    method="Gibson"
)
print("Gibson Assembly Primers:", primers)
```

**Output:**
```python
Gibson Assembly Primers: {
    'forward_primer': 'AAAAAAAAAAAAAAAAAAAAATGACCTGACTGA',
    'reverse_primer': 'TTTTTTTTTTTTTTTTTTTTTCAGGTCAT'
}
```

---

## Notes
- **Flexibility:**  
  Default settings for primer length, Tm, and homology regions can be customized by modifying the class attributes.  

- **Sequence Validation:**  
  The script ensures valid enzyme selection for Golden Gate and method selection for both approaches.

- **Applications:**  
  These primers are designed for use in experimental plasmid design pipelines specific to yeast cloning workflows.
### Yeast UTR Design: `design_utr.py`

The **UTRChooser** class is part of a Python-based pipeline designed for efficient yeast plasmid construction, focusing on the selection and optimization of untranslated regions (UTRs) for a given coding sequence (CDS). The class leverages genomic and proteomic data from *S. cerevisiae* to identify high-performing UTRs based on translation compatibility, sequence structure, and scoring criteria.

---

#### **Class Overview: UTRChooser**
The **UTRChooser** class is responsible for selecting optimal 5' and 3' UTRs based on a set of criteria that ensure efficient gene expression and compatibility with downstream molecular cloning methods. The class evaluates UTRs derived from genomic data and proteomics, and makes decisions based on edit distance, hairpin formation, and sequence checks.

---

#### **Core Features**
- **UTR Type Selection**:
  - **5' UTR**: Ensures the inclusion of a Kozak sequence for efficient ribosome assembly and translation initiation.
  - **3' UTR**: Appends a poly(A) tail to enhance mRNA stability and translation efficiency.
  
- **Scoring System**:
  - The UTRs are scored based on:
    - **Edit Distance**: Difference in the first six amino acids between the input CDS and the candidate UTR's CDS.
    - **Hairpin Count**: The number of secondary structures in the UTR that could interfere with expression.
  
- **Sequence Validation**:
  - **Forbidden Sequence Checks**: Uses a forbidden sequence checker to identify sequences that could interfere with cloning or expression (e.g., restriction sites).
  - **Kozak Sequence Compliance**: For 5' UTRs, ensures compliance with a canonical Kozak sequence, critical for translation initiation in eukaryotes.

- **Genomic and Proteomic Data Integration**:
  - Incorporates genomic data from *S. cerevisiae* and proteomics data to select high-abundance genes.
  - UTRs and CDS are derived from these genes to ensure optimal translation.

---

#### **Methods**
1. **`initiate()`**: 
   - Loads and processes genomic and proteomics data.
   - Extracts UTR and CDS sequences from high-abundance genes.
   - Creates `UTROption` instances for each potential UTR.

2. **`run(cds, end, ignores)`**:
   - Selects the best UTR for a given CDS and end type (5' or 3'), while excluding UTRs specified in the `ignores` set.
   - Returns the optimal UTR based on scoring and validation checks.

3. **`ensure_kozak()`**:
   - Validates the presence of a Kozak-like sequence in the 5' UTR for efficient translation initiation.
   - Ensures that the UTR sequence meets necessary translation initiation requirements.

4. **`forbidden_seq_check()`**:
   - Validates that the UTR does not contain forbidden sequences, such as restriction enzyme sites or other inhibitory motifs.

---

#### **Example Usage**

```python
utr_chooser = UTRChooser()
utr_chooser.initiate()

# Design a 5' UTR for the given CDS
five_prime_utr = utr_chooser.run(
    cds="ATGACCTGACTGA",  # The coding sequence
    end=5,                # 5' UTR
    ignores=set()         # No UTRs to ignore
)

# Design a 3' UTR for the given CDS
three_prime_utr = utr_chooser.run(
    cds="ATGACCTGACTGA",  # The coding sequence
    end=3,                # 3' UTR
    ignores=set()         # No UTRs to ignore
)

print("5' UTR:", five_prime_utr)
print("3' UTR:", three_prime_utr)
```

---

#### **Key Parameters**
- **`cds`**: The coding sequence for the gene for which the UTR is being selected.
- **`end`**: Specifies the UTR type, either `5` for 5' UTR or `3` for 3' UTR.
- **`ignores`**: A set of `UTROption` instances to exclude from the selection process.

---

#### **Error Handling**
- The `run()` method raises errors if:
  - The `end` parameter is not 3 or 5.
  - The CDS sequence is empty or contains invalid characters.
  - The CDS sequence length is not a multiple of 3.
  - No valid UTR options are available or remaining after filtering.

- The `ensure_kozak()` method raises an error if no Kozak-compliant 5' UTRs are found.
- The `forbidden_seq_check()` method ensures that the selected UTR does not contain any forbidden sequences that might hinder cloning or expression.

---

#### **Key Considerations**
- **Kozak Sequence**: Only valid for 5' UTRs. Ensuring the Kozak sequence is present is critical for the proper initiation of translation.
- **Forbidden Sequences**: UTRs containing certain sequences, such as restriction sites, are excluded to prevent issues during cloning.
- **Poly-A Tail**: For 3' UTRs, a poly-A tail is appended to ensure stability and efficient translation.
- **Double Pass Requirement**: The pipeline needs to be run twice: once for the 5' UTR and once for the 3' UTR, to select both UTRs for a given gene.

This pipeline streamlines the process of UTR design for yeast plasmid construction, ensuring that the selected UTRs are optimal for gene expression and compatible with downstream cloning steps.
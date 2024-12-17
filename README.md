### BioE 134 Final Project Submission: Yeast Plasmid Design

---

## **Project Overview**

This project aims to create a plasmid design pipeline for yeast, starting with inputting desired proteins and ending with experimental validation that the desired sequence can be synthesized. The pipeline includes intermediary steps such as CDS optimization, UTR and primer design, and construction file creation and validation. This project focuses on providing flexible and efficient tools for the design of yeast plasmids.

---

## **Scope of Work**

As part of the final project for BioE 134, I focused on the **UTR and Primer Design** aspects of the plasmid design pipeline. The pipeline is composed of several Python classes, each responsible for a specific step. These classes are supported by various helper functions and error-checking scripts to ensure the correctness and integrity of the design. The classes implemented are as follows:

- **UTRChooser**: A class designed to generate optimal 5' and 3' UTRs (untranslated regions) based on specific input sequences and criteria.
- **PrimerDesigner**: A class that designs primers for Golden Gate and Gibson Assembly cloning methods.
- **Helper Functions**: These include functions for reverse complementing sequences, checking for forbidden sequences in UTRs, and other tasks related to sequence validation.
  
The **error handling** mechanisms embedded in the classes ensure that incorrect or inappropriate inputs are flagged and handled, providing a robust user experience. For instance, `UTRChooser` verifies sequence lengths, ensures nucleotide characters are valid, and confirms correct 5' and 3' end specifications. Similarly, `PrimerDesigner` checks for valid restriction enzymes, cloning methods, and additional validation depending on the chosen cloning strategy.

---

## **Function Descriptions**

### **UTRChooser Class**
The `UTRChooser` class provides methods to design optimal 5' and 3' UTR sequences. It performs the following key tasks:

- **Initiate**: Initializes the default settings, including constraints for the UTR design.
- **Run**: Takes a CDS sequence and generates the corresponding UTR sequence based on the specified position (5' or 3') and additional parameters such as sequence length.

### **PrimerDesigner Class**
The `PrimerDesigner` class designs primers for Golden Gate and Gibson Assembly methods. Its functionalities include:

- **Initiate**: Initializes the primer design settings such as the target melting temperature (Tm), primer length, and restriction enzymes.
- **Run**: Designs primers based on the given CDS and UTR sequences. It accommodates both Golden Gate and Gibson Assembly methods, adjusting for overhangs and homology regions where applicable.
  
The **error handling** in this class ensures that only valid cloning methods and enzymes are used. It checks if the enzyme is part of the predefined enzyme dictionary for Golden Gate cloning, and it raises errors for unsupported methods or invalid input sequences.

### **Helper Functions**
- **Reverse Complement**: A function to return the reverse complement of a given DNA sequence.
- **Forbidden Sequence Checker**: A function to scan UTR sequences for unwanted elements such as restriction enzyme sites or secondary structure motifs.

---
## **Error Handling**

The error handling section describes the measures implemented in the pipeline to ensure robustness and prevent failures during the execution of the plasmid design process. Each class and function incorporates various checks to validate input, output, and process flow. In cases of invalid inputs or improper configurations, appropriate error messages are raised to inform the user and guide corrective actions.

### **General Error Handling**

1. **Invalid Input Sequences**: 
    - **Nucleotide Validation**: All input sequences (CDS, UTRs, primers) are checked to ensure they only contain valid nucleotide characters (`A`, `T`, `C`, `G`). If any invalid characters are detected, a `ValueError` is raised.
    - **Length Validation**: Inputs with incorrect sequence lengths (e.g., UTR length exceeding a specified maximum) will trigger a `ValueError` with a clear message indicating the issue.

2. **Method and Enzyme Validation**: 
    - For **PrimerDesigner**:
      - **Invalid Cloning Method**: The user must specify a cloning method from a predefined list (either "Golden Gate" or "Gibson"). If the method provided is not supported, a `ValueError` is raised.
      - **Unsupported Enzyme**: When using the Golden Gate method, the user must provide a valid restriction enzyme name (e.g., "BsaI"). If the enzyme is not in the predefined dictionary, an error message is raised with a list of supported enzymes.
    
3. **Sequence Compatibility**:
    - **UTR and CDS Length Compatibility**: When designing UTRs, the system checks if the length of the input sequences (5' UTR, CDS, 3' UTR) is consistent with the minimum or maximum requirements. Inconsistent sequence lengths will result in a `ValueError`.
    - **Homology Region and Overhang Sequence Lengths**: For Gibson Assembly, the homology region length is predefined. If the homology region or overhang sequences do not meet the necessary length, a `ValueError` will be raised.

4. **Invalid Region Type**:
    - **Invalid Region for Gibson Assembly**: For Gibson Assembly, primer regions can either be "upstream" or "downstream". If an invalid region type is provided, such as anything other than these two, a `ValueError` is triggered, indicating the valid region options.

### **Class-Specific Error Handling**

#### **UTRChooser Class**
The `UTRChooser` class handles errors associated with the design of UTRs:

- **Invalid 5' or 3' UTR Specification**: When generating UTRs, if the user provides an invalid specification for the position (e.g., an undefined value other than `5` or `3`), a `ValueError` is raised with a message specifying valid values.
  
- **Forbidden Sequences**: UTR sequences are scanned for forbidden sequences, such as restriction enzyme recognition sites or specific motifs. If any forbidden sequence is detected, a `ValueError` is raised with a detailed message identifying the problematic sequence.
  
- **Non-standard Sequence Characters**: UTR sequences are validated for the presence of only valid nucleotide characters (A, T, C, G). If any other characters are found, a `ValueError` is raised with a specific error message.

#### **PrimerDesigner Class**
The `PrimerDesigner` class is responsible for designing primers based on input sequences. Its error handling includes:

- **Invalid Sequence Inputs**: If the input CDS, UTR, or overhang sequences contain non-nucleotide characters or are empty, a `ValueError` is raised, specifying the invalid sequence.

- **Invalid Restriction Enzyme for Golden Gate Method**: When using Golden Gate cloning, the user must specify a valid restriction enzyme from the predefined list. If an invalid enzyme is provided, a `ValueError` is raised, listing the supported enzymes.

- **Cloning Method Validation**: The method chosen (either "Golden Gate" or "Gibson") must be valid. If an unsupported method is selected, a `ValueError` is raised, explaining the issue.

- **Inconsistent Cloning Site Sequence**: If the cloning site for Golden Gate is incompatible with the sequence or the enzyme selected, the system will raise a `ValueError`, describing the incompatibility.

- **Tm Target Calculation Errors**: During primer design, the melting temperature (Tm) is calculated for the primers. If the calculated Tm falls outside the acceptable range, the program raises a `ValueError` and advises adjusting the primer length or sequence.

- **Invalid Primer Length**: If the primer length is set to an unreasonable value (either too short or too long), the program raises an error specifying that the primer length must be within an acceptable range.

- **Invalid Region Type for Gibson Assembly**: When designing primers for Gibson Assembly, the user must specify a valid region type (`"upstream"` or `"downstream"`). If an invalid region type is provided, a `ValueError` is raised.

#### **Helper Functions Error Handling**
Several helper functions within the pipeline handle smaller tasks such as reverse complementing sequences or checking forbidden sequences. Their error handling includes:

- **Reverse Complement Validation**: If the input sequence for reverse complementing is empty or invalid (e.g., containing non-nucleotide characters), a `ValueError` will be raised.
  
- **Forbidden Sequence Checker**: When checking UTR sequences for forbidden elements, if a forbidden sequence is detected, the function will raise a `ValueError` with the sequence that triggered the error.

## **Testing**

Both the classes and their helper functions have been thoroughly tested with various input scenarios, including standard, edge, and invalid cases. A comprehensive suite of tests has been created using **pytest**, ensuring all functionality is verified and robust.

- **Test Folder**: `tests/`

The tests include:
- Valid sequences with valid UTR/Primer outputs
- Sequences containing invalid characters (non-nucleotides)
- Sequences with unsupported restriction enzymes (for primer design)
- Sequences missing valid Kozak sequences (for UTR design)
- Correct implementation of poly(A) tails (for UTR design)

This suite of tests covers a wide range of potential input errors and ensures that the tools behave as expected under various conditions.

---

## **Usage Instructions**

To use the plasmid design pipeline, clone the repository and install the required dependencies as specified in `requirements.txt`.

```bash
pip install -r requirements.txt
```

Once installed, you can use the functions from `design_primer.py` and `design_utr.py` modules as follows:

```python
from design_primer import PrimerDesigner
from design_utr import UTRChooser, UTROption

# Example CDS sequence
cds = "ATGAGTCAAGGCGGAAGAG"

# Initialize UTRChooser
utr_chooser = UTRChooser()
utr_chooser.initiate()

# Get 5' UTR
utr_5 = utr_chooser.run(cds, 5)

# Get 3' UTR while ignoring 5' UTR
utr_3 = utr_chooser.run(cds, 3, {utr_5})

# Initialize PrimerDesigner
primer_designer = PrimerDesigner()
primer_designer.initiate()

# Design Primers with Gibson Assembly
primers = primer_designer.run(cds, utr_5, utr_3, method="Gibson")
```

This will produce the primers and UTRs needed for further plasmid construction.

---

## **Conclusion**

The **UTRChooser** and **PrimerDesigner** classes form essential components of the yeast plasmid design pipeline, enabling the generation of optimized UTRs and primers necessary for efficient cloning. These tools simplify the experimental design process, ensuring that the desired sequences can be synthesized accurately for downstream applications. By incorporating error handling and input validation, these classes ensure reliable results even when given complex or non-standard input data.
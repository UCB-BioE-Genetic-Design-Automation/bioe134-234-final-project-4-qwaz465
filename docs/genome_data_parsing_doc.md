### Gene Data Extraction Documentation `genome_data_parsing.py`

This module provides functions to extract gene-related information from a GenBank file and process proteomics data. The `extract_genes_info` function parses GenBank files to retrieve gene, CDS, and UTR information, while the `proteomics_prune` function processes proteomics data to extract the most abundant genes based on a given dataset. This is to be used with S. Cerevisiae data for UTR design, however it is abstracted enough to be used with any generic GenBank/proteomics data file. Therefore, the documentation will be abstracted as such. The core of this logic has been adapted from J. Chris Anderson.

---

#### **Core Functions**

1. **`extract_genes_info(genbank_file)`**
   - **Description**: Extracts gene data from a GenBank file, including the gene name, UTR (untranslated region), and CDS (coding sequence) for each gene. It parses through all gene features in the GenBank file, and for each gene, it retrieves its corresponding CDS and UTR information.
   
   - **Parameters**:
     - `genbank_file` (str): The file path to the GenBank file to be parsed.

   - **Returns**:
     - A dictionary of gene information, indexed by the locus tag. Each entry contains:
       - `"gene"`: The gene name.
       - `"UTR"`: The 5' UTR sequence (50 bases upstream of the CDS start for forward strand genes; reverse complement for reverse strand).
       - `"CDS"`: The coding sequence for the gene.
   
   - **Example Usage**:
   ```python
   genbank_file = "data/genomic.gbff"
   genes_info = extract_genes_info(genbank_file)
   ```

---

2. **`proteomics_prune(file_path)`**
   - **Description**: Processes a proteomics data file (in tab-separated format) to extract the top 5% most abundant genes based on proteomics data. This function reads the file, converts the abundance values to numeric, handles missing values, and selects the genes with the highest abundance.
   
   - **Parameters**:
     - `file_path` (str): The file path to the proteomics data (a tab-separated file with gene abundance values).

   - **Returns**:
     - A list of tuples where each tuple contains:
       - `locus_tag`: The gene identifier.
       - `abundance`: The abundance value for that gene.

   - **Example Usage**:
   ```python
   file_path = "data/511145-WHOLE_ORGANISM-integrated.txt"
   top_5_percent_list = proteomics_prune(file_path)
   ```

---

#### **How It Works**

- **Gene Data Extraction (`extract_genes_info`)**:
  - The function uses the `SeqIO.parse` method from the `Bio` library to read and parse the GenBank file.
  - For each gene feature in the file, it extracts the gene name (`gene`), locus tag (`locus_tag`), and CDS information.
  - It then extracts the UTR sequence, which is 50 bases upstream of the CDS start for genes on the forward strand, and reverse complements the UTR for genes on the reverse strand.
  - The extracted information is stored in a dictionary for each gene, indexed by its `locus_tag`.

- **Proteomics Data Processing (`proteomics_prune`)**:
  - The function reads the proteomics data from a CSV file using the `pandas` library.
  - The `abundance` column is converted to numeric values to handle non-numeric values gracefully.
  - The data is sorted by abundance in descending order, and the top 5% of rows are selected based on the sorted abundance values.
  - The resulting list contains locus tags and their corresponding abundance values.

---

#### **Example Usage**

```python
# Extract gene information from GenBank file
genbank_file = "data/genomic.gbff"
genes_info = extract_genes_info(genbank_file)

# Process proteomics data to get top 5% abundant genes
file_path = "data/511145-WHOLE_ORGANISM-integrated.txt"
top_5_percent_list = proteomics_prune(file_path)

# Example output for the genes_info
print(genes_info["YBR098C"])  # Access UTR, CDS, and gene info for a specific gene

# Example output for the top 5% proteomics data
print(top_5_percent_list)  # Get the most abundant genes and their abundance
```

---

#### **Key Considerations**

- **GenBank File Format**: The `extract_genes_info` function assumes the GenBank file contains both `gene` and `CDS` feature types. The UTR sequence is derived based on the location of the CDS and the strand information.
  
- **Proteomics Data Format**: The `proteomics_prune` function assumes the input file has a specific tab-separated format with columns for `string_external_id` (with the format `locus_tag.xxx`) and `abundance`. The function processes and sorts the data based on the abundance values.

---

#### **Error Handling**
- **`extract_genes_info`**:
  - Raises a `ValueError` if a gene feature lacks associated CDS or UTR information.
  
- **`proteomics_prune`**:
  - Handles non-numeric `abundance` values gracefully by coercing them into `NaN`, which are dropped before further processing.
  - Raises an error if the file path is incorrect or if the file format is unexpected.

This module provides essential tools for extracting and processing genomic and proteomics data, enabling streamlined gene design and selection for further synthetic biology applications.
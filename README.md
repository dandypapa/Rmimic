# mimicR: R Interface for Peptide Mimicking

`mimicR` is an R package that provides an interface to C++ functions for generating mimic (e.g., shuffled) peptide sequences from FASTA files. This is often useful in bioinformatics, for example, in creating decoy databases for proteomics analysis.

## Installation

You can install the development version of `mimicR` from GitHub using `devtools`:

```R
# install.packages("devtools") # If you don't have devtools installed
devtools::install_github("YOUR_GITHUB_USERNAME/mimicR") # Replace with the actual path once known
```

Alternatively, you can clone the repository and build the package locally:

```bash
git clone https://github.com/YOUR_GITHUB_USERNAME/mimicR.git # Replace with actual path
R CMD build mimicR
R CMD INSTALL mimicR_<version>.tar.gz # Replace <version> with the actual version
```

## Usage

The core function of the package is `mimic_fasta()`. Here's a basic example:

```R
# Ensure the package is loaded
# library(mimicR) # Not strictly necessary if using ::

# Create a dummy FASTA file for demonstration
dummy_fasta_content <- c(">seq1_original", "PEPTIDE", ">seq2_another", "SEQUENCE")
dummy_input_path <- tempfile(fileext = ".fasta")
dummy_output_path <- tempfile(fileext = ".fasta")
writeLines(dummy_fasta_content, dummy_input_path)

# Generate mimic peptides
success <- mimicR::mimic_fasta(
  input_fasta_path = dummy_input_path,
  output_fasta_path = dummy_output_path,
  num_shuffles = 3,      # Generate 3 shuffled versions per original peptide
  replace_i = TRUE,      # Replace Isoleucine (I) with Leucine (L)
  prepend_original = TRUE # Include the original sequence in the output
)

if (success && file.exists(dummy_output_path)) {
  print(paste("Output successfully generated at:", dummy_output_path))
  # You can then read and inspect the output FASTA file
  # output_content <- readLines(dummy_output_path)
  # print(output_content)
} else {
  print("Mimic peptide generation failed or output file not found.")
}

# Clean up dummy files
if (file.exists(dummy_input_path)) file.remove(dummy_input_path)
if (file.exists(dummy_output_path)) file.remove(dummy_output_path)
```

## Key Features & Options

The `mimic_fasta()` function offers several parameters to control the peptide generation:

*   `input_fasta_path`: Path to your input FASTA file.
*   `output_fasta_path`: Path for the generated output FASTA file.
*   `min_len`: Minimum length of peptides to process.
*   `num_shuffles`: Number of mimic sequences to generate per original peptide.
*   `replace_i`: Whether to replace all Isoleucine (I) residues with Leucine (L).
*   `seed`: Seed for the random number generator for reproducible results.
*   `protein_name_prefix`: Prefix for the names of mimic proteins in the output FASTA.
*   `shared_peptide_ratio`: Ratio of shared peptides (details depend on C++ implementation).
*   `prepend_original`: If `TRUE`, the original peptide sequence is included in the output along with its mimics.
*   `infer_aa_frequency`: If `TRUE`, infers amino acid frequencies from the input file for shuffling.
*   `verbose`: Enable verbose output from the underlying C++ code.

## Vignette

For a more detailed guide and advanced usage examples, please see the package vignette:

```R
# Build vignettes if installing from source or GitHub
# devtools::install(build_vignettes = TRUE) # if you installed with devtools
# browseVignettes("mimicR")
```
*(Note: Vignette will be available after package building with vignette generation enabled.)*

## Building from Source

If you have cloned the repository, you can build the package using standard R commands:

1.  Ensure you have R development tools installed (Rtools for Windows, Xcode command-line tools for macOS, or `r-base-dev` for Linux).
2.  Navigate to the package directory in your terminal.
3.  Run `R CMD build .`
4.  Then `R CMD INSTALL mimicR_<version>.tar.gz` (replacing `<version>` with the actual version number of the built package).

To include vignettes when building from source: `R CMD build --no-manual --no-resave-data .` followed by the install command. Or ensure `VignetteBuilder: knitr` is in the `DESCRIPTION` file and build appropriately.

## License

This package is licensed under the Apache License 2.0. See the `DESCRIPTION` file for more details.
```

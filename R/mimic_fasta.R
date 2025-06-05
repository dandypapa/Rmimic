#' Generate Mimic Peptides from a FASTA File
#'
#' This function takes an input FASTA file and generates a new FASTA file
#' containing mimic (e.g., shuffled) peptides based on the original sequences.
#' It serves as an R interface to the underlying C++ 'mimic' logic.
#'
#' @param input_fasta_path Path to the input FASTA file.
#' @param output_fasta_path Path for the output FASTA file.
#' @param min_len Unsigned integer, minimum length of peptides to consider.
#'                Corresponds to the C++ '-l' option.
#' @param num_shuffles Unsigned integer, number of shuffles per original peptide.
#'                     Corresponds to the C++ '-m' (multFactor) option.
#' @param replace_i Logical, whether to replace Isoleucine (I) with Leucine (L).
#'                  Corresponds to the C++ '-I' option.
#' @param seed Unsigned integer, seed for the random number generator. 0 may indicate
#'             using a random seed. Corresponds to the C++ '-s' option.
#' @param protein_name_prefix String, prefix for protein names in the output FASTA.
#'                            Corresponds to the C++ '-p' option.
#' @param shared_peptide_ratio Double, ratio of shared peptides.
#'                             Corresponds to the C++ '-q' option.
#' @param prepend_original Logical, whether to prepend the original peptide sequence
#'                         in the output. Corresponds to the C++ '-P' option.
#' @param infer_aa_frequency Logical, whether to infer amino acid frequencies from
#'                           the input file. Corresponds to the C++ '-A' option.
#' @param verbose Logical, enable verbose output during C++ execution.
#'                Corresponds to the C++ '-v' option.
#'
#' @return Logical `TRUE` if successful and output file is created, `FALSE` otherwise.
#' @export
#' @useDynLib mimicR, .registration = TRUE
#' @importFrom Rcpp evalCpp
# This is just to ensure Rcpp is linked if sourceCpp isn't used directly.
#'
#' @examples
#' \dontrun{
#'   # Create a dummy FASTA file
#'   dummy_fasta_content <- c(">seq1", "ACGT", ">seq2", "TGCA")
#'   dummy_input_path <- tempfile(fileext = ".fasta")
#'   dummy_output_path <- tempfile(fileext = ".fasta")
#'   writeLines(dummy_fasta_content, dummy_input_path)
#'
#'   success <- mimic_fasta(
#'     input_fasta_path = dummy_input_path,
#'     output_fasta_path = dummy_output_path,
#'     num_shuffles = 2
#'   )
#'
#'   if (success && file.exists(dummy_output_path)) {
#'     print(paste("Output successfully generated at:", dummy_output_path))
#'   } else {
#'     print("Mimic peptide generation failed or output file not found.")
#'   }
#'
#'   if (file.exists(dummy_input_path)) file.remove(dummy_input_path)
#'   if (file.exists(dummy_output_path)) file.remove(dummy_output_path)
#' }
mimic_fasta <- function(input_fasta_path,
                            output_fasta_path,
                            min_len = 0,
                            num_shuffles = 1,
                            replace_i = FALSE,
                            seed = 0,
                            protein_name_prefix = "mimic|Random",
                            shared_peptide_ratio = 0.0,
                            prepend_original = FALSE,
                            infer_aa_frequency = TRUE,
                            verbose = FALSE) {

  input_fasta_path <- normalizePath(input_fasta_path, mustWork = TRUE)
  output_dir <- dirname(output_fasta_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_fasta_path <- normalizePath(output_fasta_path, mustWork = FALSE)

  success <- rcpp_mimic_fasta(
    input_fasta_path = input_fasta_path,
    output_fasta_path = output_fasta_path,
    min_len = as.integer(min_len),
    num_shuffles = as.integer(num_shuffles),
    replace_i = as.logical(replace_i),
    seed = as.integer(seed),
    protein_name_prefix = as.character(protein_name_prefix),
    shared_peptide_ratio = as.double(shared_peptide_ratio),
    prepend_original = as.logical(prepend_original),
    infer_aa_frequency = as.logical(infer_aa_frequency),
    verbose = as.logical(verbose)
  )
  return(success)
}

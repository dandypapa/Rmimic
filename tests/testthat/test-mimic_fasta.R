library(testthat)
# No need to explicitly load mimicR if tests are run by devtools::test() or R CMD check
# as it handles package loading. If running interactively, library(mimicR) would be needed.

test_that("mimic_fasta runs successfully with basic input", {
  # Create a dummy input FASTA file
  input_content <- c(">seq1", "ACDEFGHIKLMNPQRSTVWY", ">seq2", "WYVUTSRQPONMLKIHGFEDCA")
  input_file <- tempfile(fileext = ".fasta")
  output_file <- tempfile(fileext = ".fasta")
  writeLines(input_content, input_file)

  on.exit({
    if (file.exists(input_file)) file.remove(input_file)
    if (file.exists(output_file)) file.remove(output_file)
  })

  # Test with basic parameters
  success <- mimicR::mimic_fasta(
    input_fasta_path = input_file,
    output_fasta_path = output_file,
    num_shuffles = 1,
    min_len = 5 # Ensure peptides are long enough for processing
  )

  expect_true(success, info = "mimic_fasta should return TRUE on success.")
  expect_true(file.exists(output_file), info = "Output file should be created.")
  expect_gt(file.info(output_file)$size, 0, info = "Output file should not be empty.")

  # Check number of sequences in output (optional, depends on mimic logic)
  # For num_shuffles = 1, and 2 input sequences, we expect 2 output sequences if no prepend.
  # If prependOriginal = TRUE, then 4.
  # The C++ code's `multFactor_` (num_shuffles) creates that many *new* sequences.
  # Original sequences are not written unless prependOriginal_ is true.
  # Connector strings also add lines.
  # >mimic|Random_1|shuffle_1
  # PEPTIDE
  # This means 2 lines per generated peptide.
  # So, 2 input sequences * 1 shuffle/seq * 2 lines/output_peptide = 4 lines if no errors.
  # This also assumes default protein naming.

  # A simple check for lines:
  output_lines <- readLines(output_file)
  # Each peptide generates a header and a sequence line.
  # num_shuffles * number_of_input_peptides * 2 lines
  # This needs to be adjusted based on how `Peptides::run` and `cleaveProtein` work.
  # `cleaveProtein` might generate multiple peptides from one input sequence.
  # For now, let's check for a reasonable number of lines.
  # If each of the 2 input sequences (length 20) is treated as one peptide,
  # and min_len=5, then 1 shuffle each should produce 2 output peptides.
  # So, 4 lines (2 headers, 2 sequences).
  expect_gte(length(output_lines), 2 * 1 * 2,
            info = "Output file should contain expected number of lines for headers and sequences.")
})

test_that("mimic_fasta handles num_shuffles parameter", {
  input_content <- c(">seqA", "ACDEFGHIKLMNPQRSTVWY")
  input_file <- tempfile(fileext = ".fasta")
  output_file <- tempfile(fileext = ".fasta")
  writeLines(input_content, input_file)

  on.exit({
    if (file.exists(input_file)) file.remove(input_file)
    if (file.exists(output_file)) file.remove(output_file)
  })

  num_shuffles_val <- 3
  success <- mimicR::mimic_fasta(
    input_fasta_path = input_file,
    output_fasta_path = output_file,
    num_shuffles = num_shuffles_val,
    min_len = 5
  )

  expect_true(success)
  expect_true(file.exists(output_file))
  output_lines <- readLines(output_file)
  # Expected lines: 1 input seq * num_shuffles_val * 2 lines/output_peptide
  # This assumes the input sequence is treated as a single peptide to be shuffled.
  # If the input sequence "ACDEFGHIKLMNPQRSTVWY" is cleaved into multiple peptides by
  # cleaveProtein before shuffling, this count would be different.
  # The current C++ code structure in Peptides::run reads proteins, then cleaves them,
  # then shuffles the resulting peptides.
  # For simplicity in this test, we assume the input is short enough or configured
  # such that it's treated as one peptide for shuffling.
  # A more robust check would be to count ">mimic" headers.
  header_count <- sum(grepl("^>", output_lines))
  expect_equal(header_count, num_shuffles_val,
               info = "Number of output headers should match num_shuffles for a single input peptide.")
})

test_that("mimic_fasta R wrapper handles non-existent input file", {
  non_existent_file <- "this_file_does_not_exist.fasta"
  output_file <- tempfile(fileext = ".fasta")

  on.exit({
    if (file.exists(output_file)) file.remove(output_file)
  })

  # normalizePath with mustWork=TRUE in the R wrapper should catch this
  expect_error(
    mimicR::mimic_fasta(
      input_fasta_path = non_existent_file,
      output_fasta_path = output_file
    ),
    # The error message comes from normalizePath
    regexp = "cannot be found|No such file or directory"
  )
})

test_that("mimic_fasta C++ handles empty input file gracefully", {
  input_file <- tempfile(fileext = ".fasta")
  output_file <- tempfile(fileext = ".fasta")
  file.create(input_file) # Empty file

  on.exit({
    if (file.exists(input_file)) file.remove(input_file)
    if (file.exists(output_file)) file.remove(output_file)
  })

  # The C++ code might produce an empty output or an error.
  # Peptides::run() returns 0 on success.
  # An empty input might mean no peptides to process, so it could be a "success"
  # but result in an empty output file.
  # The rcpp_mimic_fasta wrapper has a check for empty output file.
  success <- mimicR::mimic_fasta(
    input_fasta_path = input_file,
    output_fasta_path = output_file
  )

  # If it's considered a success but output is empty, our wrapper returns FALSE.
  # If C++ errors out, Rcpp::stop is called.
  # If C++ returns non-zero, Rcpp::stop is called.
  # So, we expect either FALSE or an error.
  # Given the current rcpp_mimic_fasta logic, an empty output (even if run returns 0)
  # will result in FALSE (due to Rcpp::stop being called from the C++ wrapper for empty output).
  # This test needs to account for Rcpp::stop.
  expect_error(
    mimicR::mimic_fasta(
      input_fasta_path = input_file,
      output_fasta_path = output_file
    ),
    regexp = "Peptide generation reported success, but output file is missing or empty"
  )

  # The output file might be created but should be empty.
  # We need to ensure it's not created or it's empty if the function call errored out as expected.
  # If an error is thrown, the output_file might or might not exist depending on when the error occurred.
  # If it exists, it should be empty.
  if (file.exists(output_file)) {
    expect_equal(file.info(output_file)$size, 0,
                   info = "Output file, if created from empty input, should be empty.")
  }
})

# Add more tests:
# - Different min_len values
# - replace_i = TRUE
# - prepend_original = TRUE
# - Effects of seed on output (harder to test precisely without knowing exact algorithm details)
# - infer_aa_frequency = FALSE (if a mechanism to provide custom frequencies is added)

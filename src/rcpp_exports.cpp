#include <Rcpp.h>
#include "Peptides.h"
#include "AminoAcidDist.h"
#include "Option.h" // Required by Peptides.cpp
#include <fstream>
#include <set>
#include <random>
#include <vector>
#include <string>
#include <exception> // Required for std::exception

// [[Rcpp::export]]
Rcpp::LogicalVector rcpp_mimic_fasta(
    std::string input_fasta_path,
    std::string output_fasta_path,
    unsigned int min_len = 0,
    unsigned int num_shuffles = 1,
    bool replace_i = false,
    unsigned int seed = 0,
    std::string protein_name_prefix = "mimic|Random",
    double shared_peptide_ratio = 0.0,
    bool prepend_original = false,
    bool infer_aa_frequency = true,
    bool verbose = false
) {
    try {
        Peptides peptides;

        std::vector<std::string> args_str_storage;
        args_str_storage.push_back("mimicR_exec");
        args_str_storage.push_back(input_fasta_path);

        args_str_storage.push_back("-o");
        args_str_storage.push_back(output_fasta_path);

        args_str_storage.push_back("-l");
        args_str_storage.push_back(std::to_string(min_len));

        args_str_storage.push_back("-m");
        args_str_storage.push_back(std::to_string(num_shuffles));

        if (replace_i) {
            args_str_storage.push_back("-I");
        }

        args_str_storage.push_back("-s");
        args_str_storage.push_back(std::to_string(seed));

        args_str_storage.push_back("-p");
        args_str_storage.push_back(protein_name_prefix);

        args_str_storage.push_back("-q");
        args_str_storage.push_back(std::to_string(shared_peptide_ratio));

        if (prepend_original) {
            args_str_storage.push_back("-P");
        }

        if (infer_aa_frequency) {
            args_str_storage.push_back("-A");
        }

        if (verbose) {
            args_str_storage.push_back("-v");
        }

        std::vector<char*> argv_cstr;
        for (size_t i = 0; i < args_str_storage.size(); ++i) {
            argv_cstr.push_back(const_cast<char*>(args_str_storage[i].c_str()));
        }

        if (!peptides.parseOptions(static_cast<int>(argv_cstr.size()), argv_cstr.data())) {
            Rcpp::Rcerr << "Error parsing options for Peptides object in rcpp_mimic_fasta." << std::endl;
            Rcpp::stop("Failed to parse options for peptide generation.");
        }

        int result = peptides.run();

        // Check if output file was actually created and is not empty
        // This is an additional check as peptides.run() might return 0 even if output is not as expected.
        if (result == 0) {
            std::ifstream outfile_check(output_fasta_path);
            if (outfile_check.good() && outfile_check.peek() != std::ifstream::traits_type::eof()) {
                 return Rcpp::LogicalVector::create(true);
            } else {
                 Rcpp::Rcerr << "Peptides::run() reported success, but output file is missing or empty: " << output_fasta_path << std::endl;
                 Rcpp::stop("Peptide generation reported success, but output file is missing or empty: " + output_fasta_path);
            }
        } else {
             Rcpp::Rcerr << "Peptides::run() failed with result code: " << result << std::endl;
             Rcpp::stop("Peptide generation run failed with code: " + std::to_string(result));
        }

    } catch (const std::exception& e) {
        Rcpp::stop("Exception in rcpp_mimic_fasta: " + std::string(e.what()));
    } catch (...) {
        Rcpp::stop("Unknown C++ exception in rcpp_mimic_fasta.");
    }
    return Rcpp::LogicalVector::create(false); // Should not be reached
}

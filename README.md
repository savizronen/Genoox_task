# VCF File Parser

This Python script contains several classes and functions designed to process Variant Call Format (VCF) files. VCF is a standard format in bioinformatics used for storing gene sequence variations.

## Developer 

-  Ronen Saviz 

## Key Classes and Functions

- `delete_file_if_exists(filepath)`: This function deletes a file if it exists.

- `FileWriteWorker(threading.Thread)`: A worker class that uses threading to write data to a file. It's designed to improve performance when writing a large amount of data.

- `fetch_variant_details(chrom, pos, ref, alt, reference_version="hg19")`: This function fetches variant details by querying an API. If the variant details have been fetched previously, it retrieves the details from a cache instead of querying the API again.

- `classify_variant(ref, alt)`: This function classifies a variant as 'PATHOGENIC' or 'BENIGN' based on the difference between the lengths of the REF and ALT fields.

- `process_variant(sample, line, start, end, minDP, output_dir, limit, variant_counts, samples)`: This function processes a single variant from the VCF file. It fetches the variant details, classifies the variant, and writes the variant with its classification to a new VCF file.

- `VCF_parser()`: The main class of the script. It provides functionality to download a VCF file, process the file, and write the processed data to output files. It also includes methods for setting up logging, reading a configuration file, and downloading a file from an S3 bucket.

- `VCF_parser.main()`: This is the main method of the script. It downloads VCF files from specified URLs, processes the files, and writes the processed data to output files.

## Usage

This script is designed for use in a specific environment and may not work out of the box in a different environment. It expects a configuration file `VCF_parser.yaml` in a specific location and writes logs to a specific location. It also requires an API at `https://test.genoox.com/api/fetch_variant_details` to fetch variant details. 

The `VCF_parser` class uses the `process_variant` function, which is not a member of the class but is defined in the global scope of the script. This function is intended to be used as a worker function in a multithreaded environment. It processes a single variant from the VCF file, fetches the variant details, classifies the variant, and writes the variant with its classification to a new VCF file.

## Output
```
#################################### url: https://resources.genoox.com.s3.amazonaws.com/homeAssingment/demo_vcf_multisample.vcf.gz
Download the filehttps://resources.genoox.com.s3.amazonaws.com/homeAssingment/demo_vcf_multisample.vcf.gz
filename: demo_vcf_multisample.vcf.gz
File demo_vcf_multisample.vcf.gz.db not found.
############################################
start_position: None
end_position: None
min_DP: 10
limit_per_sample: 9
############################################
Open file: demo_vcf_multisample.vcf.gz and parse in parallel with workers 5  ;-)
Here is my samples: ['father', 'mother', 'proband']
Ran out of input
Processing time:  19.611694192886352539
#################################################
##################### DONE ##################### 
```
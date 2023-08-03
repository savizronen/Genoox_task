# -*- coding: utf-8 -*-
import logging.handlers
import queue
import threading
import concurrent.futures
import urllib3
from urllib.parse import urlparse
import os
import gzip
import time
from itertools import islice
import shelve
import yaml
import requests


def delete_file_if_exists(filepath):
    """
    Delete a file if it exists.
    Args:
    filepath (str): Path to the file to delete.
    Returns:
    None
    """
    try:
        os.remove(filepath)
        print(f"File {filepath} has been deleted.")
    except FileNotFoundError:
        print(f"File {filepath} not found.")


# Disable AWS certification check for public access
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


# # Open a persistent dictionary for caching
# state_db = shelve.open("cache.db", writeback=True)


class FileWriteWorker(threading.Thread):
    """
    A class that uses threading to write data to a file. This is useful for improving
    performance when writing a large amount of data.
    """

    def __init__(self, queue):
        super().__init__()
        self.queue = queue
        self.daemon = True  # Set as a daemon so it will terminate when the main thread does
        self.start()

    def run(self):
        while True:
            file_path, data = self.queue.get()
            with open(file_path, "a") as f:
                f.write(data)
            self.queue.task_done()


write_queue = queue.Queue()  # Create a queue object
write_worker = FileWriteWorker(write_queue)  # Start the file write worker

# Create a dictionary to store API responses to avoid re-querying the API for the same data
api_cache = {}


def fetch_variant_details(chrom, pos, ref, alt, reference_version="hg19"):
    """
    Fetches variant details by querying an API. If the variant details have been fetched before,
    it retrieves the details from the cache instead of querying the API again.
    """
    # Check in the cache first
    key = (chrom, pos, ref, alt, reference_version)
    if key in api_cache:
        return api_cache[key]

    # If not in the cache, query the API
    url = "https://test.genoox.com/api/fetch_variant_details"
    payload = {
        "chr": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "reference_version": reference_version
    }
    response = requests.post(url, json=payload)
    result = response.json()

    # Save the API response to the cache
    api_cache[key] = result

    return result


def classify_variant(ref, alt):
    """
    Classifies a variant as 'PATHOGENIC' or 'BENIGN' based on the difference between
    the lengths of the REF and ALT fields.
    """
    length_diff = abs(len(ref) - len(alt))
    if length_diff % 3 == 0:
        return "BENIGN"
    else:
        return "PATHOGENIC"


def process_variant(sample, line, start, end, minDP, output_dir, limit, variant_counts, samples):
    """
    Processes a single variant from the VCF file. It fetches the variant details, classifies
    the variant, and writes the variant with its classification to a new VCF file.
    """
    fields = line.strip().split("\t")
    chrom, pos, _, ref, alt, _, _, info, format_data = fields[:9]
    samples_data = fields[9:]

    format_fields = format_data.split(":")
    sample_index = samples.index(sample)
    sample_data = dict(zip(format_fields, samples_data[sample_index].split(":")))

    if "DP" in sample_data:
        try:
            dp_value = int(sample_data["DP"])
        except ValueError:
            return True
        variant_details = fetch_variant_details(chrom=chrom, pos=pos, ref=ref, alt=alt)
        gene_info = "GENE={}".format(variant_details.get("gene", ""))
        info += ";" + gene_info

        # Add classification to INFO field
        classification = classify_variant(ref, alt)
        info += ";CLASSIFICATION=" + classification

        # Create a new VCF file for the classified variant
        output_file = os.path.join(output_dir, sample, "{}_filtered.vcf".format(sample))
        data = "\t".join(fields[:7]) + "\t" + info + "\t" + "\t".join(fields[8:]) + "\n"
        write_queue.put((output_file, data))

        variant_counts[sample] += 1
        if limit is not None and variant_counts[sample] >= limit:
            return False
    return True


class VCF_parser(object):

    def __init__(self):
        self.conf = {}
        self.AppName = "VCF_parser"
        self.ConfigName = "VCF_parser.yaml"

        # Go to Home folder
        os.chdir("..")
        self.home_dir = os.getcwd() + os.sep
        self.conf_dir = self.home_dir + "conf" + os.sep
        self.tmp_dir = self.home_dir + "tmp" + os.sep
        self.logs_dir = self.home_dir + "logs" + os.sep
        self.output_dir = self.home_dir + "output" + os.sep
        # Read configuration
        self.config_file = self.conf_dir + self.ConfigName
        f = open(self.config_file, 'r')
        self.config = yaml.load(f, Loader=yaml.FullLoader)
        f.close()

        # Setup Logger
        log_file = self.logs_dir + self.AppName + ".log"
        self.log = logging.getLogger(name=self.AppName)
        self.log.setLevel(self.config["Logger"]["LogMode"])
        handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=self.config["Logger"]["maxBytes"],
                                                       backupCount=self.config["Logger"]["backupCount"],
                                                       encoding=self.config["Logger"]["encoding"])
        handler.setFormatter(logging.Formatter(self.config["Logger"]["LogFormat"]))
        self.log.addHandler(handler)

        # Print all configuration to log
        self.log.info("=================== Start ===================")
        self.log.debug("========= SETTINGS =========")
        for group in self.config.keys():
            for key in self.config[group].keys():
                self.log.debug("%s %s = %s" % (group, key, self.config[group][key]))

        os.chdir(self.tmp_dir)

        # limit (mandatory, value must be an int < 10)
        self.limit_max_default = 10
        # Number of workers
        self.max_workers = 5
        # 's3://resources.genoox.com/homeAssingment/demo_vcf_multisample.vcf.gz'
        self.urls = ['https://resources.genoox.com.s3.amazonaws.com/homeAssingment/demo_vcf_multisample.vcf.gz']

        # Refer to question
        self.start_position = self.config['Global']['start_position']
        self.end_position = self.config['Global']['end_position']
        self.min_DP = self.config['Global']['min_DP']
        self.limit_per_sample = self.config['Global']['limit_per_sample']
        self.state_db = self.config['Global']['state_db']

    def set_persistent_dictionary(self, cache_name):
        # Open a persistent dictionary for caching
        delete_file_if_exists(f"{cache_name}.db")
        self.state_db = shelve.open(f"{cache_name}.db", writeback=True)

    def parse(self, input_vcf_file):
        start_time = time.time()  # Record the start time
        self.log.info(f"start_time: {start_time}")

        # Define the input and output parameters for the script
        self.log.info(f"input_vcf_file: {input_vcf_file}")

        self.screen("############################################")
        self.screen(f"start_position: {self.start_position}")
        self.screen(f"end_position: {self.end_position}")
        self.screen(f"min_DP: {self.min_DP}")
        self.screen(f"limit_per_sample: {self.limit_per_sample}")
        self.screen("############################################")

        # Process the VCF file
        self.process_vcf_file(input_vcf_file, self.output_dir, start=self.start_position, end=self.end_position,
                              minDP=self.min_DP, limit=self.limit_per_sample)

        # Print the total processing time
        print("Processing time: ", time.time() - start_time)
        self.log.info(f"Processing time: {time.time() - start_time}")

    def screen(self, txt):
        self.log.info(txt)
        print(txt)

    def download_from_s3(self, file_url):
        """
           Download a file from an S3 bucket to a local path.
           :param file_url: The S3 HTTP URL of the file to download
       """
        try:
            self.screen(f"Download the file{file_url}")
            a = urlparse(file_url)
            filename = os.path.basename(a.path)
            response = requests.get(file_url, verify=False)
            self.screen(f"filename: {filename}")

            with open(filename, 'wb') as f:
                f.write(response.content)
            return filename
        except Exception as e:
            self.log.error("Failed to download the file, Exception:" + str(e))
            print(e)
            return None

    #  limit is mandatory
    def process_vcf_file(self, input_file, output_dir, start=None, end=None, minDP=None, limit=None):
        """
        Process the VCF file by splitting the file into samples, then processes each sample
        concurrently. It creates a new VCF file for each sample with the classified variants.
        """
        try:
            # Check mandatory
            if limit is None or limit > self.limit_max_default:
                self.log.error(f"limit - invalid value {limit}")
                return -2
            self.screen(f"Open file: {input_file} and parse in parallel with workers {self.max_workers}  ;-)")
            with gzip.open(input_file, "rt") as f, concurrent.futures.ThreadPoolExecutor(
                    max_workers=self.max_workers) as executor:
                header_lines = []
                samples = []

                # Read the VCF file header
                for line in f:
                    if line.startswith("##"):
                        header_lines.append(line)
                    elif line.startswith("#"):
                        samples = line.strip().split("\t")[9:]
                        # a list of variants (mutations) of a given genetic sample(s)
                        self.screen(f"Here is my samples: {samples}")
                        for sample in samples:
                            sample_output_dir = os.path.join(output_dir, sample)
                            os.makedirs(sample_output_dir, exist_ok=True)
                            output_file = os.path.join(sample_output_dir, "{}_filtered.vcf".format(sample))
                            with open(output_file, "w") as out_f:
                                out_f.write("".join(header_lines))
                        break

                variant_counts = {sample: 0 for sample in samples}

                start_position = self.state_db.get('position', 0)
                for line in islice(f, start_position, None):
                    for sample in samples:
                        if variant_counts[sample] < limit:
                            result = executor.submit(process_variant, sample, line, start, end, minDP, output_dir,
                                                     limit,
                                                     variant_counts, samples)
                            if not result.result():
                                samples.remove(sample)
                    # Save the current line number
                    self.state_db['position'] = start_position
                    start_position += 1
                self.state_db['position'] = start_position
            # Ensure all write tasks have completed
            write_queue.join()
            # Close the cache database
            self.state_db.close()
        except Exception as e:
            self.log.error("Failed, Exception:" + str(e))
            print(e)
            return -1

    def main(self):
        for url in self.urls:
            self.screen(f"#################################### url: {url}")
            local_file = self.download_from_s3(file_url=url)
            self.set_persistent_dictionary(local_file)
            self.parse(local_file)
            self.screen(f"#################################################")
        self.screen(f"##################### DONE #####################")


if __name__ == "__main__":
    instance = VCF_parser()
    instance.main()

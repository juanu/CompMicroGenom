#!/Users/juan/anaconda/bin/python
#Created on 7/27/2013

__author__ = 'Juan A. Ugalde'


def get_aws_credentials(credentials):
    """Read the file with the AWS credentials"""
    access_key_id = None
    secret_access_key = None

    for line in open(credentials, 'r'):
        if line.strip():
            line = line.rstrip()
            if line.startswith("User"):
                continue

            else:
                user_id, access_key_id, secret_access_key = line.split(",")

    return access_key_id, secret_access_key


###Script
if __name__ == '__main__':
    import boto
    import argparse
    import os

    program_description = "This script takes a list of protein fasta files for several genomes, makes a blastp database" \
                          "and then split the job into several spot instances on EC2. The output of each run is concatenated" \
                          "to be analyzed using orthoMCL"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-g", "--aws_credentials", type=str,
                        help="AWS credentials", required=True)
    parser.add_argument("-l","--fasta_files",type=str, help="Fasta files", required=True)
    parser.add_argument("-o", "--output_directory", type=str,
                        help="Output folder for the modified fasta files\n", required=True)

    #Add bid price as an argument in here
    bid_value = 0.3

    args = parser.parse_args()

    #Make output directory
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Get the credentials
    access_key, secret_key = get_aws_credentials(args.aws_credentials)

    #AMI = ami-1ad03273

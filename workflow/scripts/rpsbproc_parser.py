import os
import csv
from tqdm import tqdm
import pandas as pd

import logging
import colored_logging

import defaults

def parse_rpsbproc_output(input_dir, output_dir):
    domain_hits = []
    fieldnames = [
        'File_Name',  # Added the File_Name field as the first column
        'Session_ordinal', 'Program', 'Version', 'Database', 'Score_matrix', 'Evalue_threshold',
        'Query_ID', 'Seq_type', 'Seq_length', 'Definition',
        'Hit_type', 'PSSM_ID', 'From', 'To', 'Evalue', 'Bit_score',
        'Accession', 'Short_name', 'Incomplete', 'Superfamily_PSSM_ID'
    ]

    # Get list of files to process
    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]
    total_files = len(files)
    file_bar = tqdm(total=total_files, desc='Processing files...')

    for filename in files:
        input_file = os.path.join(input_dir, filename)
        if not os.path.isfile(input_file):
            file_bar.update(1)
            continue  # Skip directories

        logging.info(f"Processing file: {input_file}...")
        data_started = False
        current_session = {}
        current_query = {}

        # Get the file name without extension
        file_name_without_ext = os.path.splitext(filename)[0]

        with open(input_file, 'r') as file:
            lines = file.readlines()

        index = 0
        total_lines = len(lines)

        while index < total_lines:
            line = lines[index].strip()

            # Skip header lines until DATA
            if not data_started:
                if line == 'DATA':
                    data_started = True
                index += 1
                continue

            # Process lines after DATA
            if line.startswith('SESSION'):
                # Start of a new session
                parts = line.split('\t')
                current_session = {
                    'Session_ordinal': parts[1],
                    'Program': parts[2],
                    'Version': parts[3],
                    'Database': parts[4],
                    'Score_matrix': parts[5],
                    'Evalue_threshold': parts[6]
                }
                index += 1

            elif line.startswith('QUERY'):
                # Start of a new query
                parts = line.split('\t')
                current_query = {
                    'Query_ID': parts[1],
                    'Seq_type': parts[2],
                    'Seq_length': parts[3],
                    'Definition': parts[4]
                }
                index += 1

                # Check if the next line is DOMAINS
                if index < total_lines and lines[index].strip() == 'DOMAINS':
                    index += 1  # Skip 'DOMAINS' line
                    # Process domain hits until ENDDOMAINS
                    while index < total_lines and not lines[index].strip().startswith('ENDDOMAINS'):
                        domain_line = lines[index].strip()
                        domain_parts = domain_line.split('\t')
                        domain_hit = {
                            'File_Name': file_name_without_ext,  # Add the file name to the domain hit
                            'Session_ordinal': domain_parts[0],
                            'Query_ID': domain_parts[1],
                            'Hit_type': domain_parts[2],
                            'PSSM_ID': domain_parts[3],
                            'From': domain_parts[4],
                            'To': domain_parts[5],
                            'Evalue': domain_parts[6],
                            'Bit_score': domain_parts[7],
                            'Accession': domain_parts[8],
                            'Short_name': domain_parts[9],
                            'Incomplete': domain_parts[10],
                            'Superfamily_PSSM_ID': domain_parts[11],
                            'Seq_type': current_query['Seq_type'],
                            'Seq_length': current_query['Seq_length'],
                            'Definition': current_query['Definition'],
                            'Program': current_session['Program'],
                            'Version': current_session['Version'],
                            'Database': current_session['Database'],
                            'Score_matrix': current_session['Score_matrix'],
                            'Evalue_threshold': current_session['Evalue_threshold']
                        }
                        domain_hits.append(domain_hit)
                        index += 1
                    index += 1  # Skip 'ENDDOMAINS'
                else:
                    # No domain hits for this query
                    no_domain_hit = {
                        'File_Name': file_name_without_ext,  # Add the file name to the no domain hit
                        'Session_ordinal': current_session.get('Session_ordinal', ''),
                        'Program': current_session.get('Program', ''),
                        'Version': current_session.get('Version', ''),
                        'Database': current_session.get('Database', ''),
                        'Score_matrix': current_session.get('Score_matrix', ''),
                        'Evalue_threshold': current_session.get('Evalue_threshold', ''),
                        'Query_ID': current_query.get('Query_ID', ''),
                        'Seq_type': current_query.get('Seq_type', ''),
                        'Seq_length': current_query.get('Seq_length', ''),
                        'Definition': current_query.get('Definition', ''),
                        'Hit_type': '',
                        'PSSM_ID': '',
                        'From': '',
                        'To': '',
                        'Evalue': '',
                        'Bit_score': '',
                        'Accession': '',
                        'Short_name': '',
                        'Incomplete': '',
                        'Superfamily_PSSM_ID': ''
                    }
                    domain_hits.append(no_domain_hit)
                # Skip 'ENDQUERY' if present
                if index < total_lines and lines[index].strip() == 'ENDQUERY':
                    index += 1
            else:
                index += 1

        file_bar.update(1)

    file_bar.close()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write to CSV with progress bar
    csv_output_file = os.path.join(output_dir, 'domains.csv')
    logging.info(f"Writing CSV to {csv_output_file}")
    total_hits = len(domain_hits)
    with open(csv_output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        hits_bar = tqdm(total=total_hits, desc='Writing CSV', leave=False)
        for hit in domain_hits:
            writer.writerow(hit)
            hits_bar.update(1)
        hits_bar.close()
    logging.info(f"CSV output written to {csv_output_file}")

    # Write to Parquet
    logging.info("Writing Parquet file")
    df = pd.DataFrame(domain_hits)
    parquet_output_file = os.path.join(output_dir, 'domains.parquet')
    df.to_parquet(parquet_output_file, index=False)
    logging.info(f"Parquet output written to {parquet_output_file}")

if __name__ == '__main__':
    colored_logging.colored_logging(log_file_name='rpsbproc_parser.txt')

    parse_rpsbproc_output(defaults.RPSBPROC_OUTPUT_DIR, defaults.TABLE_OUTPUT_DIR)

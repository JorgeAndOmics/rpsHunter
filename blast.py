from tqdm import tqdm
import subprocess
import pandas as pd
import os

import defaults

import logging
from colored_logging import colored_logging

def blaster(command: str, input_database_path, query_file_path, subject, outfmt: str):
    """
    Runs a BLAST search for a given object against a given database

        Parameters
        ----------

            :param command: The command to run BLAST.
            :param input_database_path: The path to the database.
            :param query_file_path: The path to the query file.
            :param subject: The particular genome against whose database it's being BLASTed
            :param outfmt: The output format for the BLAST results. Default is 5.

        Returns
        -------
            :returns: The output of the BLAST search, captured from std_out.

        Raises
        ------
            :raise Exception: If an error occurs while running BLAST.
    """
    input_path = os.path.join(input_database_path, subject, subject) if subject else input_database_path
    subject = subject or input_database_path
    try:
        # Construct the BLAST command
        blast_command = [
            command,
            '-db',
            os.path.join(input_path),
            '-query',
            query_file_path,  # BLAST+ doesn't take in SeqRecord objects, but files
            '-evalue',
            str(defaults.E_VALUE),
            '-outfmt',
            outfmt,
        ]

        # Run the command
        # logging.debug(f"Running {command} for {instance.accession}: {instance.probe} probe against '{subject}'"
        #               f"\n{instance.display_info()}")
        result = subprocess.run(blast_command, capture_output=True, text=True)

        # Output captured in result.stdout
        blast_output = result.stdout
        if blast_output.strip():
            return blast_output

        logging.warning('BLAST output is empty')
        return None
    except Exception as e:
        logging.error(f'An error occurred while running {command}: {str(e)}')
        return None

def parse_blast_output(blast_output):
    """Parses the BLAST output from string into a list of lists for DataFrame construction."""
    # Ensure the output is not empty
    if not blast_output:
        return []

    # Split the output into rows by newline
    rows = blast_output.strip().split('\n')

    return [row.split('\t') for row in rows]


if __name__ == '__main__':
    # Set up logging
    colored_logging(log_file_name='blast.txt')

    # Define the column headers for the BLAST output
    columns = ['Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
               'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score', 'Species']

    df = pd.DataFrame(columns=columns)

    with tqdm(total=len(defaults.SPECIES)) as pbar:
        for species in defaults.SPECIES:
            pbar.set_description(f'Processing {species.replace("_", " ")}...')
            # Run the BLAST command and capture output
            blast_output = blaster(command='tblastn',
                                   input_database_path=defaults.SPECIES_DB,
                                   query_file_path=defaults.QUERY_SEQ,
                                   subject=species,
                                   outfmt='6')

            if parsed_blast_output := parse_blast_output(blast_output=blast_output):
                # Append the species to the parsed BLAST output
                for row in parsed_blast_output:
                    row.append(species)

                # Append the parsed BLAST output to the DataFrame
                new_df = pd.DataFrame(parsed_blast_output, columns=columns)
                df = pd.concat([df, new_df], axis=0)
            else:
                logging.warning('No data returned from BLAST.')

            pbar.update(1)

    # Save the DataFrame to PARQUET and CSV file
    output_csv_path = os.path.join(defaults.TABLE_OUTPUT_DIR, 'blast.csv')
    output_parquet_path = os.path.join(defaults.TABLE_OUTPUT_DIR, 'blast.parquet')

    df.to_csv(output_csv_path, index=False)
    logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

    df.to_parquet(output_parquet_path, index=False)
    logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')
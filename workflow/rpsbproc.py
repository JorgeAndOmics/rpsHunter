import os
import subprocess

import defaults


def main(input_dir, output_dir, db_path, t_option):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # List all .asn files in the input directory
    input_files = [f for f in os.listdir(input_dir) if f.endswith('.asn')]

    for filename in input_files:
        input_file = os.path.join(input_dir, filename)
        output_filename = f'{os.path.splitext(filename)[0]}.txt'
        output_file = os.path.join(output_dir, output_filename)

        cmd = [
            'rpsbproc',
            '-i', input_file,
            '-d', db_path,
            '-t', t_option,
            '-o', output_file
        ]

        print(f'Processing {input_file} -> {output_file}')
        subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main(input_dir=defaults.ASN_RPSBLAST_DIR,
         output_dir=defaults.RPSBPROC_OUTPUT_DIR,
         db_path=defaults.RPSBPROC_DB,
         t_option='doms')

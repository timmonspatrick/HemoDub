# HemoDub

## Usage
1. Run scrape_databases.py with Python 3 to download the HTML files of the websites
2. Run read_files.py with Python 3 to read the HTML files and produce peptides.json
3. read_files.py relies heavily on comment_parser.py to parse comments describing hemolytic activities.
4. Run desciptors.py with Python 2 to fully describe the peptide sequences contained in peptides.json and produce the numpy array peptide_array.npy

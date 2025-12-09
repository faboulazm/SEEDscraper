# SEED SCRAPER

A webscraper pipeline to automate extracting proteins seed sequences from UniProtKB database using query entries. Stores concatenated sequences into fasta file and SQL database. 

## Dependencies: 
1. requests
2. biopython
3. beautifulsoup4

## Installation & Use: 
1. `git clone https://github.com/faboulazm/SEEDscraper.git`
2. `cd SEEDscraper`
3. `python src/SEEDscraper.py`

## Workflow: 
1. Run the pipeline & enter queries for proteins of interest. 
2. Enter EXIT to terminate query entry. 
3. Concatenated FASTA & protein.db will be created. 

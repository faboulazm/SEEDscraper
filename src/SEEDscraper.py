from bs4 import BeautifulSoup
from urllib.parse import quote
from Bio import SeqIO
import requests
import sqlite3

all_fasta_data = ""
print("=====================================")
print("Welcome to Fatima's SEEDscraper :)")
print("Enter UniProt queries one by one. Type EXIT to finish and create database.")
print("Query example: MurA AND (taxonomy_name:2) AND (reviewed:true)")
print("=====================================")

while True: 
    query = input("Enter your query & type EXIT to finish: ")
    
    if query.strip().upper() == "EXIT":
        break

    encoded_query = quote(query)
    url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded_query}&format=fasta"
    response = requests.get(url)

    if response.status_code == 200:
        fasta_data = response.text
        print(f"Fetched {len(fasta_data.splitlines()) // 2} sequences (approx)")
        all_fasta_data += fasta_data + "\n" 
    else: 
        print(f"Failed to retrieve data. Status code: {response.status_code }")
print("=====================================")
print("Finished retrieving sequences. ")

with open("all_sequences.fasta", "w") as f:
          f.write(all_fasta_data)
print("======================================")
print("Sequences saved to all_sequences.fasta")
print("======================================")

from io import StringIO
fasta_io = StringIO(all_fasta_data)


sequences = []

from io import StringIO

fasta_io = StringIO(all_fasta_data)

for record in SeqIO.parse(fasta_io, "fasta"):
    uniprot_id = record.id.split("|")[-1]  # P12345
    description = record.description
    sequence = str(record.seq)
    
    
    parts = description.split(" ")
    protein_name = parts[1] if len(parts) > 1 else ""
    organism = ""
    for part in parts:
        if part.startswith("OS="):
            organism = part.replace("OS=", "")
    
    sequences.append((uniprot_id, protein_name, organism, sequence))

print(f"Parsed {len(sequences)} sequences.")


conn = sqlite3.connect("proteins.db")
c = conn.cursor()

c.execute("""
CREATE TABLE IF NOT EXISTS proteins (
    uniprot_id TEXT PRIMARY KEY,
    protein_name TEXT,
    organism TEXT,
    sequence TEXT
)
""")
conn.commit()


for seq in sequences:
    try:
        c.execute(
            "INSERT INTO proteins (uniprot_id, protein_name, organism, sequence) VALUES (?, ?, ?, ?)",
            seq
        )
    except sqlite3.IntegrityError:
        pass  

conn.commit()
print("=====================================")
print(f"{len(sequences)} sequences stored in proteins.db")
print("=====================================")
conn.close()


from Bio import SeqIO

gb = SeqIO.parse(open("bs168genbank.gbff", "r"), "genbank")

for record in gb:
    for items in record:
        print(record.features)

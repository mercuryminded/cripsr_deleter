#making a script that will find genes in a genome of a specific annotation and then design CRISPR guide RNAs to delete
#those genes

'''genbank and fasta file downloaded from here https://www.ncbi.nlm.nih.gov/genome/665?genome_assembly_id=300274'''

from Bio import SeqIO

gb_file = SeqIO.parse(open("bs168genbank.gbff", "r"), "genbank")

counter = 0
for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
    feats = gb_record.features
    print("Name: {}, Features: {}".format(gb_record.name, len(feats)))
    for item in feats:
        print(item.gene)
        counter += 1
        if counter > 20:
            break
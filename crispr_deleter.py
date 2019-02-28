#making a script that will find genes in a genome of a specific annotation and then design CRISPR guide RNAs to delete
#those genes

'''genbank and fasta file downloaded from here https://www.ncbi.nlm.nih.gov/genome/665?genome_assembly_id=300274'''

from Bio import SeqIO

gb_file = SeqIO.parse(open("bs168genbank.gbff", "r"), "genbank")

counter = 0

#trouble with using qualifiers is that the tags are different for all the entries. Need to have exceptions for when tag
#is not there, or the error will stop the code at the first entry.

for record in gb_file:
    print(repr(record.seq))
    for i in range(1,10):
        feature = record.features[i]
        #print(type(feature))
        #print(feature.qualifiers)
        for key in feature.qualifiers['gene']:
            print(key)



'''
for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
    feats = gb_record.features
    print("Name: {}, Features: {}".format(gb_record.name, len(feats)))
    for item in feats:
        print(item.gene)
        counter += 1
        if counter > 20:
            break
'''
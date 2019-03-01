#making a script that will find genes in a genome of a specific annotation and then design CRISPR guide RNAs to delete
#those genes

'''genbank and fasta file downloaded from here https://www.ncbi.nlm.nih.gov/genome/665?genome_assembly_id=300274'''

from Bio import SeqIO

gb = SeqIO.parse(open("bs168genbank.gbff", "r"), "genbank")

#trouble with using qualifiers is that the tags are different for all the entries. Need to have exceptions for when tag
#is not there, or the error will stop the code at the first entry.

#find feature.qualifiers['function'] contains 'Biosynthesis' and feature.qualifiers['product'] contains 'antibiotic'
#or something like that

def finds_tags(genbank_file):
    feature_list=[]
    gb_file = SeqIO.parse(open(genbank_file, "r"), "genbank")
    for record in gb_file:
        print(record.name)
        all_feats = record.features
        for feature in all_feats:
            y = 0
            for qualifier in feature.qualifiers:
                if qualifier == 'gene':
                    y += 1
                if qualifier == 'product':
                    y += 1
                if qualifier == 'function':
                    y += 1
            if y == 3:
                if 'antibiotic' in feature.qualifiers['product'][0]:
                    l = feature.location
                    feature_list.append([feature.qualifiers['product'][0], int(l.start), int(l.end), l.strand])
                    print(feature.qualifiers['product'][0])
                    print(l.start)
                    print(l.end)
                    print(l.strand)
    print(feature_list)
                #need to find a way to save the specific feature in an easy to access way

finds_tags("bs168genbank.gbff")
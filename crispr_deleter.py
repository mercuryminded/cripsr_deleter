# making a script that will find genes in a genome of a specific annotation and then design CRISPR guide RNAs to delete
# those genes

'''genbank and fasta file downloaded from here https://www.ncbi.nlm.nih.gov/genome/665?genome_assembly_id=300274'''

from Bio import SeqIO
import pickle

gb = SeqIO.parse(open("bs168genbank.gbff", "r"), "genbank")

# trouble with using qualifiers is that the tags are different for all the entries. Need to have exceptions for when tag
# is not there, or the error will stop the code at the first entry.
# find feature.qualifiers['function'] contains 'Biosynthesis' and feature.qualifiers['product'] contains 'antibiotic'
# or something like that

# make a thing that shows all of the possible search terms for functions
def lists_terms(genbank_file, search_field, search_term, setting):
    """settings: 0 = list, 1 = hunt"""
    gb_file = SeqIO.parse(open(genbank_file, "r"), "genbank")
    file_name = str(genbank_file[0:4]) + '_' + str(search_term) + '_' + str(search_field)
    file_object = open(file_name, 'wb')
    for record in gb_file:
        feature_set = set()
        feature_locations = []
        for feature in record.features:
            y = 0
            for qualifier in feature.qualifiers:
                # checks that the feature is a gene, has a listed product and listed function
                if qualifier == 'gene':
                    y += 1
                if qualifier == 'product':
                    y += 1
                if qualifier == 'function':
                    y += 1
                    # print(feature.qualifiers['function'])
            if y == 3 and setting == 0:
                for entry in feature.qualifiers[search_field]:
                    feature_set.add(entry)
            if y == 3 and setting == 1:
                for entry in feature.qualifiers[search_field]:
                    if search_term in entry:
                        record.features.extract(record.seq)
                        feature_locations.append([feature.qualifiers[search_field], int(l.start), int(l.end), l.strand])
    if setting == 0:
        print("Here is the set of features")
        print(feature_set)
    if setting == 1:
        pickle.dump(feature_locations, file_object)
        file_object.close()


def collects_terms():
    """insert docstring here"""
    x = ['function', 'product']
    file_name = input("Input name of genbank file (make sure it's in the working directory): ")
    search_field_no = 5
    while search_field_no != 0 and search_field_no != 1:
        try:
            search_field_no = int(input('Search in: 0 = function, 1 = product: '))
        except ValueError:
            print("Do it right lmao")
    if search_field_no == 0:
        lists_terms(file_name, x[search_field_no],'no_term', 0)
        search_term = input('Type in one of the functions listed above: ')
    else:
        search_term = input('Name of gene product to search for: ')
    # possible modification: allow multiple search terms
    lists_terms(file_name, x[search_field_no], search_term, 1)

collects_terms()

"""
def finds_tags(term_tup):
    feature_list = []
    genbank_file = term_tup[0]
    search_term = term_tup[1]
    search_field = term_tup[2]
    gb_file = SeqIO.parse(open(genbank_file, "r"), "genbank")

    file_name = str(genbank_file[0:4]) + '_' + str(search_term) + '_' + str(search_field)
    file_object = open(file_name, 'wb')

    for record in gb_file:
        print(record.name)
        for feature in record.features:
            y = 0
            for qualifier in feature.qualifiers:
                # checks that the feature is a gene, has a listed product and listed function
                if qualifier == 'gene':
                    y += 1
                if qualifier == 'product':
                    y += 1
                if qualifier == 'function':
                    y += 1
                    # print(feature.qualifiers['function'])
            if y == 3:
                for entry in feature.qualifiers[search_field]:
                    if search_term in entry:
                        l = feature.location
                        feature_list.append([feature.qualifiers[search_field], int(l.start), int(l.end), l.strand])
    pickle.dump(feature_list, file_object)
    file_object.close()
    # makes a list [name, start, end, strand]
    # need to find a way to save the specific feature in an easy to access way

# finds_tags("bs168genbank.gbff", 'Biosynthesis', 'function')

finds_tags(collects_terms())
"""

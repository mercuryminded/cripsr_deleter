# making a script that will find genes in a genome of a specific annotation and then design CRISPR guide RNAs to delete
# those genes

'''genbank and fasta file downloaded from here https://www.ncbi.nlm.nih.gov/genome/665?genome_assembly_id=300274'''

from Bio import SeqIO
import pickle
import os

gb = SeqIO.parse(open("bs168genbank.gbff", "r"), "genbank")

# trouble with using qualifiers is that the tags are different for all the entries. Need to have exceptions for when tag
# is not there, or the error will stop the code at the first entry.
# find feature.qualifiers['function'] contains 'Biosynthesis' and feature.qualifiers['product'] contains 'antibiotic'
# or something like that

# make a thing that shows all of the possible search terms for functions


def lists_terms(genbank_file, search_field, search_term, setting):
    """settings: 0 = list, 1 = search"""
    gb_file = SeqIO.read(open(genbank_file, "r"), "genbank")
    pickle_name = str(genbank_file[0:3]) + '_pickle'
    feature_set = set()
    feature_dict = {}
    if setting == 0 and os.path.isfile(pickle_name) and search_field == 'function':
        unpick = pickle.load(open(pickle_name, 'rb'))
        print("Here is the set of features")
        print(*unpick, sep="\n")
    for feature in gb_file.features:
        y = 0
        for qualifier in feature.qualifiers:
            # checks that the feature is a gene, has a listed product and listed function
            if qualifier == 'gene':
                y += 1
            elif qualifier == 'product':
                y += 1
            elif qualifier == 'function':
                y += 1
                # print(feature.qualifiers['function'])
        if setting == 0 and not os.path.isfile(pickle_name):
            if y == 3 and setting == 0: # a feature must have above 3 properties to be searched
                for entry in feature.qualifiers[search_field]:
                    feature_set.add(entry)
        if y == 3 and setting == 1:
            for entry in feature.qualifiers[search_field]:
                if search_term in entry:
                    feature_dict[feature.qualifiers['product'][0]] = feature.extract(gb_file.seq)
    if setting == 0 and not os.path.isfile(pickle_name):
        print("Here is the set of features")
        print(feature_set, sep="\n")
        file_object = open(pickle_name, 'wb')
        pickle.dump(feature_set, file_object)
    if setting == 1:
        f = open(genbank_file[0:3] + '_' + search_field + search_term, "w+")
        f.write(feature_dict)
        print(feature_dict)


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


#lists_terms('bs168genbank.gbff', 'function', 'no_term', 0)
collects_terms()

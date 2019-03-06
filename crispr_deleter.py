# making a script that will find genes in a genome of a specific annotation and then design CRISPR guide RNAs to delete
# those genes

'''genbank and fasta file downloaded from here https://www.ncbi.nlm.nih.gov/genome/665?genome_assembly_id=300274'''

from Bio import SeqIO
import pickle
import os

# gb = SeqIO.parse(open("bs168genbank.gbff", "r"), "genbank")
# trouble with using qualifiers is that the tags are different for all the entries. Need to have exceptions for when tag
# is not there, or the error will stop the code at the first entry.
# find feature.qualifiers['function'] contains 'Biosynthesis' and feature.qualifiers['product'] contains 'antibiotic'
# or something like that

# make a thing that shows all of the possible search terms for functions


def lists_terms(genbank_file, search_field, search_term, setting):
    """settings: 0 = list, 1 = search"""
    gb_file = SeqIO.read(open(genbank_file, "r"), "genbank")
    pickled_set_name = str(genbank_file[0:3]) + '_set_pickle'
    pickled_dict_name = str(genbank_file[0:3]) + '_' + search_term + '_dict'
    feature_set = set()
    feature_dict = {}
    if setting == 0 and os.path.isfile(pickled_set_name) and search_field == 'function':
        unpick_set = pickle.load(open(pickled_set_name, 'rb'))
        print("Here is the set of features")
        print(*unpick_set, sep="\n")
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
        if setting == 0 and not os.path.isfile(pickled_set_name):
            if y == 3 and setting == 0: # a feature must have above 3 properties to be searched
                for entry in feature.qualifiers[search_field]:
                    feature_set.add(entry)
        if y == 3 and setting == 1:
            for entry in feature.qualifiers[search_field]:
                if search_term in entry:
                    feature_dict[feature.qualifiers['product'][0]] = feature.extract(gb_file.seq)
    if setting == 0 and not os.path.isfile(pickled_set_name):
        print("Here is the set of features")
        print(feature_set, sep="\n")
        file_object = open(pickled_set_name, 'wb')
        pickle.dump(feature_set, file_object)
        file_object.close()
    if setting == 1:
        file_object1 = open(pickled_dict_name, 'wb')
        pickle.dump(feature_dict, file_object1)
        file_object1.close()
        # print(feature_dict)
        return feature_dict


def collects_terms():
    """insert docstring here"""
    x = ['function', 'product']
    file_name = input("Input name of genbank file (make sure it's in the working directory): ")
    search_field_no = 5
    while search_field_no != 0 and search_field_no != 1:
        try:
            search_field_no = int(input('Type 0 to search by function (list provided), type 1 to search by product name: '))
        except ValueError:
            print("Type in 0 or 1 to choose: ")
    if search_field_no == 0:
        lists_terms(file_name, x[search_field_no],'no_term', 0)
        search_term = input('Type in one of the functions listed above: ')
    else:
        search_term = input('Name of gene product to search for: ')
    # possible modification: allow multiple search terms
    return lists_terms(file_name, x[search_field_no], search_term, 1)


# lists_terms('bs168genbank.gbff', 'function', 'no_term', 0)
# collects_terms()


def startup_sequence():
    """docstring"""
    setting = int(input("Type 0 for new input, type 1 for opening saved file: "))
    if setting == 0:
        seq_dict = collects_terms()
        makes_rna(seq_dict)
    elif setting == 1:
        x = input("Type in name of pickled file: ")
        seq_dict = pickle.load(open(x, 'rb'))
        makes_rna(seq_dict)

def makes_rna(dict_feats):
    """docstring"""
    position_dict = {}
    for name, sequence in dict_feats.items():
        position_list = []
        for x in range(len(sequence)):
            if sequence[x+1:x+3] == 'GG':
                position_list.append(x)
            if sequence.complement()[x+1:x+3] == 'GG':
                position_list.append(x*(-1))
            position_dict[name] = position_list
    for name, position in position_dict.items():
        if len(position) > 0:
            print(name + '\n' + str(position))
        else:
            print(name + '\nThere is no available PAM site')

startup_sequence()
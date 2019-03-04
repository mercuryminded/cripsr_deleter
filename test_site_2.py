from Bio import SeqIO
import pickle
gb = SeqIO.read(open("bs168genbank.gbff", "r"), "genbank")

for item in gb.features:
    print(item)
    for x in item.qualifiers:
        print(type(x))


print(len(gb.features))
print(type(gb.seq))


def lists_terms(genbank_file, search_field, search_term, setting):
    """settings: 0 = list, 1 = hunt"""
    gb_file = SeqIO.read(open(genbank_file, "r"), "genbank")
    file_name = str(genbank_file[0:4]) + '_' + str(search_term) + '_' + str(search_field)
    file_object = open(file_name, 'wb')
    feature_set = set()
    feature_locations = []
    print(1)
    for feature in gb_file.features:
        y = 0
        print(2)
        for qualifier in feature.qualifiers:
            print(2.1)
            # checks that the feature is a gene, has a listed product and listed function
            if qualifier == 'gene':
                y += 1
            elif qualifier == 'product':
                y += 1
            elif qualifier == 'function':
                y += 1
                # print(feature.qualifiers['function'])
        if y == 3 and setting == 0:
            print(3)
            for entry in feature.qualifiers[search_field]:
                feature_set.add(entry)
        if y == 3 and setting == 1:
            for entry in feature.qualifiers[search_field]:
                if search_term in entry:
                    feature.extract(gb_file.seq)
                    feature_locations.append([feature.qualifiers[search_field], int(l.start), int(l.end), l.strand])
    if setting == 0:
        print("Here is the set of features")
        print(feature_set)
    if setting == 1:
        pickle.dump(feature_locations, file_object)
        file_object.close()

lists_terms('bs168genbank.gbff', 'function', 'no_term', 0)

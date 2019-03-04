import pickle
from Bio import SeqIO
gb_file = SeqIO.parse(open('bs168genbank.gbff', "r"), "genbank")


def makes_rna(pickles):
    """unpacks into lists with the """
    sequence_list = []
    file_object = open(pickles, 'rb')
    x = pickle.load(file_object)
    for record in gb_file:
        for item in x:
            print(item)
            start = item[1]
            stop = item[2]

    print(sequence_list)


makes_rna("bs16_Circulate_function")

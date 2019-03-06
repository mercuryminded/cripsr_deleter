import pickle
from Bio import SeqIO
gb_file = SeqIO.parse(open('bs168genbank.gbff', "r"), "genbank")


x = 'bs1_Biosynthesis_dict'

seq_dict = pickle.load(open(x, 'rb'))
for key, value in seq_dict.items():
    print(value)
    print(value.complement())
    print(value.complement().transcribe() + 'UUU')


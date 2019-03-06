import pickle
from Bio import SeqIO
gb_file = SeqIO.parse(open('bs168genbank.gbff', "r"), "genbank")


x = 'bs1_Biosynthesis_dict'

search_term = 'tang'

search_term = [x.strip() for x in search_term.split(',')]

print(search_term)
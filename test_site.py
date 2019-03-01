
import pickle
def makes_rna(pickles):
    file_object = open(pickles, 'r')
    x = pickle.load(file_object)
    print(x)

makes_rna("bs16_helicase_product")
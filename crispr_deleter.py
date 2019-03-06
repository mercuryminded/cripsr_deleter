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


def startup_sequence():
    """Boots up the script and will either load a saved feature dictionary or makes a new set of sequences.

    Takes no arguments, requires user input:
        First input is the setting: 0 or 1
            0: Calls collects_terms() to make a new set of sequences
            1: Loads a file (ending in _dict, containing a dictionary of {names:coding sequences})
                Second user input is the file name for setting 1
        Both settings will call finds_pam_sites() and then makes_rna()
    """
    setting = int(input("Type 0 for new input, type 1 for opening saved file: "))
    # both blocks below do the same thing, except
    # setting 0 calls collects_terms() and gets new sequences
    # setting 1 retrieves a saved dictionary of sequences
    if setting == 0:
        seq_dict = collects_terms()
        dict_positions = finds_pam_sites(seq_dict) # calls finds_pam_sites() to make dictionary of name:positions
        makes_rna(dict_positions, seq_dict) # uses the dictionary of sequences and dictionary of positions
    elif setting == 1:
        x = input("Type in name of pickled file: ")
        if not os.path.isfile('./' + x):
            raise ValueError
        seq_dict = pickle.load(open(x, 'rb'))
        dict_positions = finds_pam_sites(seq_dict)
        makes_rna(dict_positions, seq_dict)


def collects_terms():
    """Asks user for input regarding file name and search terms
    Takes no arguments, requires user input:
        First input is the file name of a genbank file in the current working directory
        Second input is

    """
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
    elif search_field_no == 1:
        lists_terms(file_name, x[search_field_no], 'no_term', 0)
        search_term = input('List of possible products above, type keywords separated by \',\' to search: ')
    # possible modification: allow multiple search terms
    if "," in search_term:
        search_term = [x.strip() for x in search_term.split(',')]
        # input will now accept multiple search terms and make them into lists
    return lists_terms(file_name, x[search_field_no], search_term, 1)


def lists_terms(genbank_file, search_field, search_term, setting):
    """settings: 0 = list, 1 = search"""
    gb_file = SeqIO.read(open(genbank_file, "r"), "genbank")
    pickled_set_name = str(genbank_file[0:3]) +'_'+search_field+ '_set_pickle'
    pickled_dict_name = str(genbank_file[0:3]) + '_' + search_term[0].lower() + str(len(search_term)) + '_dict'
    feature_set = set()
    feature_dict = {}
    if setting == 0 and os.path.isfile(pickled_set_name):
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
                search_score = 0
                for x in search_term:
                    if x.lower() in entry.lower():
                        search_score += 1
                        if search_score == len(search_term):
                            feature_dict[feature.qualifiers['product'][0]] = feature.extract(gb_file.seq)
                    # this extracts the sequence from the specific feature and adds it as a value with product name as key
    if setting == 0 and not os.path.isfile(pickled_set_name):
        print("Here is the set of features")
        print(*feature_set, sep="\n")
        file_object = open(pickled_set_name, 'wb')
        pickle.dump(feature_set, file_object)
        file_object.close()
    if setting == 1:
        file_object1 = open(pickled_dict_name, 'wb')
        pickle.dump(feature_dict, file_object1)
        file_object1.close()
        # print(feature_dict)
        return feature_dict


# lists_terms('bs168genbank.gbff', 'function', 'no_term', 0)
# collects_terms()


def finds_pam_sites(dict_feats):
    """Finds target cut sites (PAM sequences) for SpCas9. Checks on both strands and
    outputs a negative number for the non-coding strand."""
    position_dict = {}
    for name, sequence in dict_feats.items():
        position_list = []
        for x in range(len(sequence)):
            if sequence[x+1:x+3] == 'GG':
                position_list.append(x)
            if sequence.reverse_complement()[x+1:x+3] == 'GG':
                # reverse complement find GG. Will output negative number
                position_list.append(x*(-1))
            position_dict[name] = position_list
    for name, pos in position_dict.items():
        if len(pos) == 0:
            position_dict[name] = 'No PAM site available'
    return position_dict


def makes_rna(dict_positions, dict_feats):
    """docstring"""
    while not (0 <= promoter_setting <= 3):
        try:
            promoter_setting = int(input('0: Bacillus subtilis promoter, 1: Escherichia coli promoter, 2: Mammalian promoter, 3: No promoter'))
        except ValueError:
            print("Type one of the numbers listed to choose: ")
    rna_dict = {}
    unavailable_dict = {}
    scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT'
    promoter_dict = {
        0: 'AAA',
        1: 'CCC',
        2: 'GGG',
        3: 'TTT'
    }
    for name, position in dict_positions.items():
        if type(position) == str:
            print('Removing ' + name + ' from dictionary')
            unavailable_dict[name] = dict_positions[name] # adds unavailable proteins to a dictionary
            del dict_positions[name] # removes them from the dictionary
    for name, position in dict_positions.items():
        seq = dict_feats[name]
        while name not in rna_dict.keys():
            for entry in position:
                if 20 <= entry <= 50:
                    rna_dict[name] = (promoter_dict[promoter_setting] + seq[entry-20:entry] + scaffold)
                    break
                elif -len(seq)+20 <= entry <= -len(seq)+50:
                    rna_dict[name] = (promoter_dict[promoter_setting] + seq[entry-20:entry] + scaffold)
                    break
                elif 51 <= entry <= len(seq)/2:
                    rna_dict[name] = (promoter_dict[promoter_setting] + seq[entry-20:entry] + scaffold)
                    break
                elif -len(seq)+51 <= entry <= -len(seq)/2:
                    rna_dict[name] = (promoter_dict[promoter_setting] + seq[entry-20:entry] + scaffold)
                    break
                elif position.index(entry) == len(position)-1 and not rna_dict[name]:
                    unavailable_dict[name] = 'No PAM site available within first half of CDS'
    for key, value in rna_dict.items():
        print(key)
        print(value.complement())
        print('\n')
    print("These sequences are now ready to be cloned into an appropriate DNA vector for use.")
    if unavailable_dict:
        for key, value in unavailable_dict.items():
            print(key)
            print(value)
            print('\n')


startup_sequence()
"""
A script that will find genes in a genome with the specified annotation and then generate CRISPR guide RNA sequences to
target those genes for breakage by SpCas9. NHEJ then repairs the double strand break and has a chance of causing
a frame shift error. The target sites (PAM sites) are picked to be as near to the start of the coding sequence so that the
frame shift error will cause a significant portion of the coding sequence to be altered.
"""

# genbank and fasta file downloaded from here https://www.ncbi.nlm.nih.gov/genome/665?genome_assembly_id=300274

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pickle
import os


def startup_sequence():
    """
    Boots up the script and will either load a saved feature dictionary or makes a new set of sequences.

    Takes no arguments, requires user input:
        First input is the setting: 0 or 1
            0: Calls collects_terms() to make a new set of sequences
            1: Loads a file (ending in _dict, containing a dictionary of {names:coding sequences})
                Second user input is the file name for setting 1
        Both settings will call finds_pam_sites() and then makes_rna()
    """
    setting = int(input("Type 0 for new input, type 1 for opening a saved feature dictionary: "))
    # both blocks below do the same thing, except
    # setting 0 calls collects_terms() and gets new sequences
    # setting 1 retrieves a saved dictionary of sequences
    if setting == 0:
        seq_dict = collects_terms()  # calls collects_terms() to make a dictionary of {'name':sequences}
        positions_dict = finds_pam_sites(seq_dict)  # calls finds_pam_sites() to make dictionary of {'name':[positions]}
        makes_rna(positions_dict, seq_dict)  # uses the dictionary of sequences and dictionary of positions
    elif setting == 1:
        x = ''
        while not os.path.isfile(x):
            try:
                x = input("Type in name of pickled file (ends with '_pickle'): ")
                seq_dict = pickle.load(open(x, 'rb'))
            except FileNotFoundError:
                print('Type in the name of a genbank file in the current working directory (including file extension).')
        positions_dict = finds_pam_sites(seq_dict)
        makes_rna(positions_dict, seq_dict)


def collects_terms():
    """
    Asks user for input regarding file name and search terms
    Takes no arguments, requires user input:
        First input is the file name of a genbank file in the current working directory
        Second input is choosing to search in functions of genes (0) or their products (1)
            List of all possible search terms for the chosen field will be printed for reference
    :return lists_terms(): returns the result of lists_terms() on setting 1, which returns a dictionary of
                           {'Name':Sequence} pairs
    """
    x = ['function', 'product']
    file_name = ''
    while not os.path.isfile(file_name):
        try:
            file_name = input("Input name of genbank file (make sure it's in the working directory): ")
            gb_file = SeqIO.read(open(file_name, "r"), "genbank")
            print("Genbank file ID is: " + gb_file.id)
        except FileNotFoundError:
            print("Please make sure the genbank file is in the same directory as the python script, "
                  "and that you have typed in the full name of the file")
    search_field_no = 5
    while search_field_no != 0 and search_field_no != 1:
        try:
            search_field_no = int(input('Type 0 to search by function, '
                                        'type 1 to search by product name: '))
            search_field = x[search_field_no]
        except ValueError:
            print("Type in 0 or 1 to choose: ")
        except IndexError:
            print("Type in 0 or 1 to choose: ")
    if search_field == 'function':
        lists_terms(file_name, search_field,'no_term', 0)
        search_term = input('List of possible products above, type in one or more keywords to search. '
                            "Separate multiple keywords with a ','. Please note that coding sequences usually have "
                            "only one or two functions: ")
    elif search_field == 'product':
        lists_terms(file_name, search_field, 'no_term', 0)
        search_term = input('List of possible products above, type in one or more keywords to search. '
                            "Separate multiple keywords with a ',' \n : ")
    search_term = [x.strip() for x in search_term.split(',')]
    # input will now accept multiple search terms and make them into lists
    return lists_terms(file_name, search_field, search_term, 1)


def lists_terms(genbank_file, search_field, search_term, setting):
    """
    Takes four arguments returned by collects_terms(), pickles a set of all functions and products if not already done.
    Searches the genome for coding sequences annotated with the search terms and adds them into feature_dict.
    Then pickles a feature_dict file to load from. Finally returns feature_dict for finds_pam_sites() to use.
    :param genbank_file: String. Name of the genbank file containing a genome to be searched
    :param search_field:
    :param search_term:
    :param setting:
    :return feature_dict:
    """
    gb_file = SeqIO.read(open(genbank_file, "r"), "genbank")
    pickled_set_name = str(genbank_file[0:3]) + '_' + search_field + '_set_pickle'
    pickled_feature_dict = str(genbank_file[0:3]) + '_' + search_term[0] + '_' + str(len(search_term)) + '_dict'
    feature_set = set()
    feature_dict = {}
    if setting == 0 and os.path.isfile(pickled_set_name):  # checks if the file is present
        unpick_set = pickle.load(open(pickled_set_name, 'rb'))  # unpickles and loads the set
        print(*unpick_set, sep="\n")
        print("Above is the set of features to choose from. Please wait... ")
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
        if setting == 0 and not os.path.isfile(pickled_set_name):
            # setting 0 generates initial reference list of elements
            if y == 3 and setting == 0:  # a feature must have above 3 properties to be searched
                for entry in feature.qualifiers[search_field]:
                    feature_set.add(entry)
        if y == 3 and setting == 1:
            # setting 1 generates the dictionary of sequences of the search results
            for entry in feature.qualifiers[search_field]:
                search_score = 0
                for x in search_term:
                    if x.lower() in entry.lower():
                        search_score += 1
                        if search_score == len(search_term):
                            feature_dict[feature.qualifiers['product'][0]] = feature.extract(gb_file.seq)
                            # extracts the sequence from the feature and adds it as a value with product name as key
    if setting == 0 and not os.path.isfile(pickled_set_name):
        print(*feature_set, sep="\n")
        print("Above is the set of features to choose from. Please wait... ")
        file_object = open(pickled_set_name, 'wb')
        pickle.dump(feature_set, file_object)
        file_object.close()
    if setting == 1:
        file_object1 = open(pickled_feature_dict, 'wb')
        pickle.dump(feature_dict, file_object1)
        return feature_dict


def finds_pam_sites(dict_feats):
    """
    Finds target cut sites (PAM sequences) for SpCas9. Checks on both strands and
    outputs a negative number for the non-coding strand.
    :param dict_feats: Dictionary of {'Name':Sequence} pairs
    :return position_dict: Dictionary of {'Name':PAM site position} pairs
    """
    position_dict = {}
    for name, sequence in dict_feats.items():
        position_list = []
        position_list_negative = []
        for x in range(len(sequence)):
            if sequence[x+1:x+3] == 'GG' and x >= 20:
                position_list.append(x)
            if sequence.complement()[x+1:x+3] == 'GG' and x >= 20:
                # reverse complement find GG on complementary strand
                position_list_negative.append(-x)
        while name not in position_dict.keys():
            try:
                reverse = int(max(position_list_negative))  # -25
                forward = int(min(position_list))  # 28
                if abs(reverse) < forward:
                    position_dict[name] = int(reverse)
                else:
                    position_dict[name] = int(forward)
            except ValueError:  # occurs when one or both of the position_lists has no values
                # removes the list if empty and continues with one, or marks as unavailable if bth are empty
                if len(position_list) == 0 and len(position_list_negative) == 0:
                    position_dict[name] = 'No PAM site available'
                elif len(position_list) > 0 and len(position_list_negative) == 0:
                    position_dict[name] = int(min(position_list))
                elif len(position_list) == 0 and len(position_list_negative) > 0:
                    position_dict[name] = int(max(position_list_negative))
                else:
                    position_dict[name] = 'No PAM site available/Unknown error occurred'
    return position_dict


def makes_rna(dict_positions, dict_feats):
    """
    This function generates a FASTA file containing DNA sequences encoding the guide RNAs in FASTA format.
    If any sequences do not have a suitable site for targeting, the names of the sequences will be printed.
    Takes one user input which is the name of the FASTA file to be saved.
    :param dict_positions: Dictionary containing {'Name':[positions]} pairs
    :param dict_feats: Dictionary containing {'Name':Sequence} pairs
    :return:
    """
    rna_dict = {}
    unavailable_dict = {}
    pam_dict = {}
    scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT'
    for name, position in dict_positions.items():
        if type(position) == str:
            print('Removing ' + name + ' from dictionary')
            unavailable_dict[name] = dict_positions[name] # adds unavailable proteins to a dictionary
            del dict_positions[name] # removes them from the dictionary
    for name, position in dict_positions.items():
        seq = dict_feats[name]
        while name not in rna_dict.keys():
            if position < 0 and abs(position) <= len(seq)/2:
                # checks if position is on negative strand and uses complementary sequence if true
                rna_dict[name] = seq.complement()[abs(position)-20:abs(position)] + scaffold
                pam_dict[name] = position-1
                break
            elif 0 < position <= len(seq)/2:  # only considers the position if it is in the first half of the CDS
                rna_dict[name] = (seq[position-20:position] + scaffold)  # collects the gRNA target sequence
                pam_dict[name] = position-1  # marks the position of the PAM site to be used
                break
            else:
                unavailable_dict[name] = 'No PAM site available within first half of CDS'
    save_name = input("Name of output fasta file (do not include file extension): ")
    with open(save_name + ".fasta", "w+") as output_handle:
        for key, value in rna_dict.items():
            value_r = SeqRecord(value)
            # converts value into a SeqRecord file to be handled by SeqIO
            value_r.id = value_r.description = 'CRISPR SpCas9 gRNA targeting ' + key + ' at position ' \
                                               + str(pam_dict[key]) + ' of the coding sequence'
            SeqIO.write(value_r, output_handle, 'fasta')
        output_handle.close()
    print("These sequences are now ready to be cloned into an appropriate DNA vector with a promoter for use.")
    if len(unavailable_dict) > 0:
        for key, value in unavailable_dict.items():
            print(key)
            print(value)
            print('\n')
        print("item(s) listed above do not have suitable PAM sequences for SpCas9 knockout.")

startup_sequence()
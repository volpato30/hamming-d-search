NodeCount = 0


class TrieNode:
    def __init__(self):
        self.peptide = None
        self.children = {}

        global NodeCount
        NodeCount += 1

    def insert(self, peptide: list):
        node = self
        for amino_acid in peptide:
            if amino_acid not in node.children:
                node.children[amino_acid] = TrieNode()

            node = node.children[amino_acid]

        node.peptide = peptide


def build_trie_from_list(peptide_list: list):
    root_trie = TrieNode()
    for peptide in peptide_list:
        root_trie.insert(peptide)
    return root_trie


def build_trie_dict_from_file(file_name: str) -> dict:
    trie_dict = {}
    peptide_dict = {}
    with open(file_name, 'r') as f:
        for line in f:
            peptide = list(line.strip('\n'))
            if not peptide:
                # skip empty lines
                continue
            length = len(peptide)
            if length not in peptide_dict:
                peptide_dict[length] = [peptide]
            else:
                peptide_dict[length].append(peptide)
    for length, peptide_list in peptide_dict.items():
        trie_dict[length] = build_trie_from_list(peptide_list)
    return trie_dict


def build_trie_dict_from_denovo_file(file_name: str) -> dict:
    import csv
    trie_dict = {}
    peptide_dict = {}
    with open(file_name, 'r') as fr:
        reader = csv.reader(fr, delimiter='\t')
        header = next(reader)
        seq_index = header.index("predicted_sequence")
        for line in reader:
            if not line[seq_index]:
                continue
            peptide = line[seq_index].split(',')
            peptide = [x[0] for x in peptide]  # drop mod
            length = len(peptide)
            peptide = ''.join(peptide)
            if length not in peptide_dict:
                peptide_dict[length] = [peptide]
            else:
                peptide_dict[length].append(peptide)
    for length, peptide_list in peptide_dict.items():
        trie_dict[length] = build_trie_from_list(peptide_list)
    return trie_dict


def aa_match(aa1, aa2):
    return aa1 == aa2


def aa_match_ignore_I_L(aa1, aa2):
    if aa1 == "I" and aa2 == "L":
        return True
    elif aa1 == "L" and aa2 == "I":
        return True
    else:
        return aa1 == aa2


def search_hamming_d(root_trie: TrieNode, peptide: list, d: int, aa_match_function=aa_match):
    result = []
    current_index = 0
    current_cost = 0
    for amino_acid in root_trie.children:
        _search_recursive(root_trie.children[amino_acid], amino_acid, peptide,
                          current_index, current_cost, d, result, aa_match_function)
    return result


def _search_recursive(node: TrieNode, amino_acid, peptide: list, current_index: int,
                      current_cost: int, d: int, result: list, aa_match_function):
    if not aa_match_function(amino_acid, peptide[current_index]):
        current_cost += 1
    if current_cost > d:
        return
    current_index += 1

    # stopping criteria
    if node.peptide:
        result.append(node.peptide)
        return

    for aa in node.children:
        _search_recursive(node.children[aa], aa, peptide, current_index, current_cost, d, result, aa_match_function)


def same_peptide(pep1:list, pep2:list):
    assert len(pep1) == len(pep2)
    for i, a in enumerate(pep1):
        if not aa_match_ignore_I_L(a, pep2[i]):
            return False
    return True


if __name__ == '__main__':
    ## search wild type peptides
    trie_dict = build_trie_dict_from_file('./identified.peptides')
    output_file = 'wild_type_matched.peptides'
    with open(output_file, 'w') as fw:
        with open('./denovo.peptides', 'r') as fr:
            for i, line in enumerate(fr):
                peptide = list(line.strip())
                length = len(peptide)
                if len(peptide) < 8 or len(peptide) > 14:
                    continue
                result_1 = search_hamming_d(trie_dict[length], peptide=peptide, d=1,
                                          aa_match_function=aa_match_ignore_I_L)
                if result_1:
                    flag = True
                    for peptide_ in result_1:
                        # if a peptide has an exact match in the idenfied peptide list,
                        # then do not include it in neoantigen list
                        if same_peptide(peptide, peptide_):
                            flag = False
                    if flag:
                        fw.write(''.join(peptide) + '\n')
                if i % 1000 == 0:
                    print(f"processing {i}th peptide")


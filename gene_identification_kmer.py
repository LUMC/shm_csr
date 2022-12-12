#!/usr/bin/env python3

from collections import Counter
from typing import Dict, Iterator, List, Tuple


def kmer_profile(sequence: str, k: int):
    kmers = (sequence[i:i+k] for i in range(len(sequence) - k))
    return Counter(kmers)


def immuno_gene_profiles():
    ca = "catccccgaccagccccaaggtcttcccgctgagcctctgcagcacccagccagatgggaacgtggtcatcgcctgcctgg"
    cg = "ctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctctgggggcacagcggcc"
    ce = "gcctccacacagagcccatccgtcttccccttgacccgctgctgcaaaaacattccctcc"
    cm = "gggagtgcatccgccccaacc"
    # lambda/kappa referesearchstringsnce sequence variable nucleotides
    ca1_mutations = {38: 't', 39: 'g', 48: 'a', 49: 'g', 51: 'c', 68: 'a',
                     73: 'c'}
    ca2_mutations = {38: 'g', 39: 'a', 48: 'c', 49: 'c', 51: 'a', 68: 'g',
                     73: 'a'}
    cg1_mutations = {0: 'c', 33: 'a', 38: 'c', 44: 'a', 54: 't', 56: 'g',
                     58: 'g', 66: 'g', 132: 'c'}
    cg2_mutations = {0: 'c', 33: 'g', 38: 'g', 44: 'g', 54: 'c', 56: 'a',
                     58: 'a', 66: 'g', 132: 't'}
    cg3_mutations = {0: 't', 33: 'g', 38: 'g', 44: 'g', 54: 't', 56: 'g',
                     58: 'g', 66: 'g', 132: 'c'}
    cg4_mutations = {0: 't', 33: 'g', 38: 'g', 44: 'g', 54: 'c', 56: 'a',
                     58: 'a', 66: 'c', 132: 'c'}

    def mutate(sequence: str, mutations: Dict[int, str]):
        return "".join(c if i not in mutations else mutations[i]
                       for i, c in enumerate(sequence))

    iga1 = mutate(ca, ca1_mutations)


def generate_sequence_and_id_from_summary(summary_file: str
                                          ) -> Iterator[Tuple[str, str]]:
    with open(summary_file, "rt") as summary:
        header = next(summary)
        column_names = header.strip("\n").split("\t")
        id_column = column_names.index("Sequence ID")
        sequence_column = column_names.index("Sequence")
        for line in summary:
            values = line.strip("\n").split("\t")
            id = values[id_column]
            try:
                sequence = values[sequence_column]
            except IndexError:  # weird rows without a sequence
                sequence = ""
            yield id, sequence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",
                        help="The 1_Summary file from an IMGT zip file")
    parser.add_argument("--output",
                        help="The annotated output file to be merged back "
                             "with the summary file")
    args = parser.parse_args()

    with open(args.output, "wt") as output:
        output.write("Sequence ID\tbest_match\tnt_hit_percentage\t"
                     "chunk_hit_percentage\tstart_locations\n")
        for id, sequence in generate_sequence_and_id_from_summary(args.input):
            best_match, subclass_hits, class_hits, start_locations = \
                match_sequence(sequence, compiledregex)
            variable_nucs = subclass_vars[best_match]
            if variable_nucs:
                subclass_percentage = round(subclass_hits * 100 /
                                            variable_nucs)
            else:
                subclass_percentage = 100
            class_percentage = round(class_hits * 100 / class_chunks[best_match])
            output.write(f"{id}\t{best_match}\t{subclass_percentage}\t"
                         f"{class_percentage}\t{start_locations}\n")

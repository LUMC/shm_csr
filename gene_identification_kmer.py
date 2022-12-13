#!/usr/bin/env python3

import argparse
from typing import Dict, Iterator, List, Set, Tuple


def generate_kmers(sequence: str, k: int) -> Iterator[str]:
    return (sequence[i:i+k] for i in range(len(sequence) - k))


def distinct_kmers(sequence: str, k: int) -> Set[str]:
    return set(generate_kmers(sequence, k))


def immuno_gene_profiles(k: int) -> Tuple[Dict[str, Set[str]],
                                          Dict[str, Dict[str, Set[str]]]]:
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

    ca1 = mutate(ca, ca1_mutations)
    ca2 = mutate(ca, ca2_mutations)
    cg1 = mutate(cg, cg1_mutations)
    cg2 = mutate(cg, cg2_mutations)
    cg3 = mutate(cg, cg3_mutations)
    cg4 = mutate(cg, cg4_mutations)

    iga1_profile = distinct_kmers(ca1, k)
    iga2_profile = distinct_kmers(ca2, k)
    igg1_profile = distinct_kmers(cg1, k)
    igg2_profile = distinct_kmers(cg2, k)
    igg3_profile = distinct_kmers(cg3, k)
    igg4_profile = distinct_kmers(cg4, k)
    ige_profile = distinct_kmers(ce, k)
    igm_profile = distinct_kmers(cm, k)
    # Combine all subclass profiles to provide the best matching for the
    # entire class
    iga_profile = iga1_profile | iga2_profile
    igg_profile = igg1_profile | igg2_profile | igg3_profile | igg4_profile
    class_profiles = {
        "IGA": iga_profile,
        "IGE": ige_profile,
        "IGG": igg_profile,
        "IGM": igm_profile,
    }
    subclass_profiles = {
        "IGA": {
            "IGA1": iga1_profile,
            "IGA2": iga2_profile,
        },
        "IGG": {
            "IGG1": igg1_profile,
            "IGG2": igg2_profile,
            "IGG3": igg3_profile,
            "IGG4": igg4_profile,
        }
    }
    return class_profiles, subclass_profiles


def match_sequence(sequence, class_profiles: Dict[str, Set[str]],
                   subclass_profiles: Dict[str, Dict[str, Set[str]]]):
    # retrieve k from the length of the first item in the set for IGA
    k = len(next(iter(class_profiles["IGA"])))
    profile = distinct_kmers(sequence, k)
    matches = [
        (len(profile & gene_profile), gene)
        for gene, gene_profile in class_profiles.items()
    ]
    # Sort so the highest match is on top
    matches.sort(reverse=True)
    match_score, best_gene = matches[0]
    if match_score == 0:
        return "IGA1", 0, 0
    subclasses = subclass_profiles.get(best_gene)
    if not subclasses:
        return best_gene, 1, match_score / len(class_profiles[best_gene])
    subclass_matches = [
        (len(profile & subclass_profile), subclass)
        for subclass, subclass_profile in subclasses.items()
    ]
    subclass_matches.sort(reverse=True)
    subclass_match, best_subclass = subclass_matches[0]
    max_matches = len(subclasses[best_subclass])
    subclass_score = subclass_match / max_matches
    # Use min so there is not above 100% score
    class_score = min(match_score, max_matches) / max_matches
    return best_subclass, class_score, subclass_score


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
    parser.add_argument("-k", "--kmer-size", help="The size of k",
                        type=int, default=7)
    args = parser.parse_args()

    class_profiles, subclass_profiles = immuno_gene_profiles(args.kmer_size)

    with open(args.output, "wt") as output:
        output.write("Sequence ID\tbest_match\tnt_hit_percentage\t"
                     "chunk_hit_percentage\tstart_locations\n")
        for id, sequence in generate_sequence_and_id_from_summary(args.input):
            best_match, subclass_score, class_score = \
                match_sequence(sequence, class_profiles, subclass_profiles)
            output.write(f"{id}\t{best_match}\t{round(subclass_score * 100)}\t"
                         f"{round(class_score * 100)}\t[]\n")

if __name__ ==  "__main__":
    main()

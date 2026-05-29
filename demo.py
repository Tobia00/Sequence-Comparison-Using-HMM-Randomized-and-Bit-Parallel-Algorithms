import os
from randomized_alignment import read_fasta, demo as randomized_demo
from pair_hmm_alignment import demo as pair_hmm_demo
from bit_parallel_alignment import demo as bit_parallel_demo


def main():
    base_dir = os.path.dirname(__file__)
    seq1_path = os.path.join(base_dir, "data", "e_coli_demo.fasta")
    seq2_path = os.path.join(base_dir, "data", "salmonella_demo.fasta")

    seq1 = read_fasta(seq1_path)
    seq2 = read_fasta(seq2_path)

    print("Biological Sequence Alignment Demo")
    print("=" * 60)
    print("Sequence 1 length:", len(seq1))
    print("Sequence 2 length:", len(seq2))

    randomized_demo(seq1, seq2)
    pair_hmm_demo(seq1, seq2)
    bit_parallel_demo(seq1[:14], seq2[:14])


if __name__ == "__main__":
    main()

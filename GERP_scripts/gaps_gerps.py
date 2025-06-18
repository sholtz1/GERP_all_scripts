from Bio import SeqIO
from collections import deque
from dataclasses import dataclass
import sys
import argparse

 
def parse_args(args):
    # Instantiate the parser
    parser = argparse.ArgumentParser(description="Make gerp file gapped")

    # Required positional argument
    parser.add_argument("--alignment_path", type=str, help="Path to alignment file")
    parser.add_argument("--reference_name", type=str, help="Name of reference sequence in alignment")
    parser.add_argument("--gerp_file", type=str, help="Path to gerp file")
    parser.add_argument("--gerp_output", type=str, help="Path to output gerp file")
    parser.add_argument("--chromosome", type=int, help="Chromosome number to print to output file")
    return parser.parse_args()


@dataclass
class GerpStuff:
    alignment_path: str
    reference_name: str
    gerp_file: str
    gerp_output: str
    chromosome: int
    
    def __post_init__(self):
        self.get_ipython_alignment()
        self.get_gap_idx()
        self.get_gerp()
        self.write_gerp_gapped()

    def get_ipython_alignment(self):
        self.alignment = SeqIO.to_dict(SeqIO.parse(self.alignment_path, "fasta"))
    
    def get_gap_idx(self):
        self.gap_idx = [1 if base == "N" else 0 for idx, base in enumerate(self.alignment[self.reference_name].seq)]

    def get_gerp(self):
        gerp_scores = []
        with open(self.gerp_file, "r") as f:
            for line in f:
                if line:
                    gerp_scores.append(line.strip())
        self.gerp_scores = deque(gerp_scores)

    def write_gerp_gapped(self):
        assert len(self.gap_idx) == len(self.alignment[self.reference_name].seq)
        gerp_gap = []
        for idx, score in enumerate(self.gap_idx):
            if score == 1:
                gerp_gap.append("NA\tNA")
            else:
                x = self.gerp_scores.popleft()
                gerp_gap.append(x)
        with open(f"{self.gerp_output}", "w") as f:
            for idx, score_str in enumerate(gerp_gap):
                print(f"{self.chromosome}\t{idx}\t{score_str}", file=f)    


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    GerpStuff(**vars(args))

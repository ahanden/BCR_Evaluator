#!/usr/bin/python
"""
Bisulfite conversion ratio estimator.

BCR evaluator is used to evaluate bisulfite conversion ratio in WGSBS
experiment.

Authors:
    Quanyuan He Ph.D (hqyone@hunnu.edu.com)

Contributors:
    Adam Handen (adam.handen@gmail.com)

Licensed under the MIT License
"""

from argparse import ArgumentParser
from re import compile as re_compile

from pysam import FastxFile


class BCReval():
    """Bisulfite conversion rate calculator."""

    c_blocks = [
        'TTTTAA',
        'CTTTAA',
        'TCTTAA',
        'TTCTAA',
        'CCTTAA',
        'TCCTAA',
        'CTCTAA',
        'CCCTAA',
    ]
    g_blocks = [
        'TTAAAA',
        'TTAAAG',
        'TTAAGA',
        'TTAGAA',
        'TTAAGG',
        'TTAGGA',
        'TTAGAG',
        'TTAGGG',
    ]
    patterns = [
        re_compile('(([CT]{3}TAA)+)'),  # Patterns of C strand
        re_compile('((TTA[GA]{3})+)'),  # Patterns of G strand
    ]

    def __init__(self, strand, min_rep):
        """
        Create a new BCReval object.

        Parameters:
            strand (int): Expected sequence strandedness
            min_rep (int): Minimum telomeric repeats to record
        """
        self.strand = strand
        self.frequencies = {}
        self.min_rep = min_rep

    def sum_frequencies(self, blocks):
        """
        Count and weight the processed telomeric sequences.

        Parameters:
            blocks (list): Telomeric patterns (keys for frequencies)

        Returns:
            (int) Weighted sum of telomeric sequences
        """
        weights = [
            (1, 1),
            (2, 1),
            (3, 1),
            (4, 2),
            (5, 1),
            (6, 3),
            (7, 3),
        ]
        freq = self.frequencies
        freq_gen = (freq.get(blocks[blk], 0) * wt for blk, wt in weights)
        return sum(freq_gen)

    def bcr(self):
        """
        Compute the bisulfite conversion rate.

        Returns:
            (float) BCR
        """
        if self.strand is None:
            return (self._bcr(0) + self._bcr(1)) / 2
        return self._bcr(self.strand)

    def longest_block(self, seq, pattern):
        """
        Find the longest telomeric repeat.

        Parameters:
            seq (str): DNA sequence
            pattern (re.pattern): Telomeric pattern to search for

        Returns:
            (str) Longest telomeric repeat
        """
        matches = pattern.findall(seq)
        if not matches:
            return None
        return max([_[0] for _ in matches], key=len)

    def process_seq(self, seq, strand=None):
        """
        Process an individual DNA sequence.

        Parameters:
            seq (str): DNA sequence
            strand (int): Seq strandedness
        """
        if strand is None:
            self.process_seq(seq, 0)
            self.process_seq(seq, 1)
            return
        block = self.longest_block(seq, self.patterns[strand])
        if block is None:
            return
        if len(block) / 6 < self.min_rep:
            return
        for _ in range(0, len(block), 6):
            sb = block[_:_ + 6]
            self.frequencies[sb] = self.frequencies.get(sb, 0) + 1

    def proc_fastq(self, fastq_path):
        """
        Process a FastQ file.

        Parameters:
            fastq_path (str): Path to FastQ formatted file
        """
        with FastxFile(fastq_path) as fq:
            for read in fq:
                self.process_seq(read.sequence, self.strand)

    def _bcr(self, strand):
        telos = self.c_blocks if self.strand == 0 else self.g_blocks
        total = sum(self.frequencies.get(_, 0) for _ in telos)
        if total == 0:
            return 0
        return self.sum_frequencies(telos) / (3 * total)


def main(fastq_path, strand=None, min_length=5):
    """
    Primary method for commandline tool.

    Parameters:
        fastq_path (str): Path to FastQ file
        strand (int): Strandedness of the fastq file (default None)
        min_length (int): Minimum number of telomeric repeats
    """
    counter = BCReval(strand, min_length)
    counter.proc_fastq(fastq_path)
    print(counter.bcr())


if __name__ == '__main__':
    parser = ArgumentParser(
        prog='BCReval.py',
        description='Estimates bisulfite conversion rate from telomere reads.',
    )
    parser.add_argument('filename')
    parser.add_argument(
        '-s',
        '--strand',
        choices=[0, 1],
        type=int,
        help='Whether to run in stranded mode (options are 0 and 1)',
    )
    parser.add_argument(
        '-n',
        '--min-length',
        default=5,
        type=int,
        help='Minimum telomeric repeats (default is 5)',
    )
    args = parser.parse_args()
    main(args.filename, args.strand, args.min_length)

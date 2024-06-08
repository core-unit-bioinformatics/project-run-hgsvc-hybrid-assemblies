#!/usr/bin/env python3

import argparse as argp
import collections as col
import enum
import functools as fnt
import math
import pathlib as pl
import re
import sys

import pandas as pd


class CigarOp(enum.Enum):
    IDENTITY = 0
    MISMATCH = 1
    MATCH = 2
    INSERTION = 3
    DELETION = 4
    HARDMASK = 5
    SOFTMASK = 6


class CigarNormalizer:

    def __init__(self):
        self.known_ops = {
            "=": CigarOp.IDENTITY,
            "X": CigarOp.MISMATCH,
            "M": CigarOp.MATCH,
            "I": CigarOp.INSERTION,
            "D": CigarOp.DELETION,
            "H": CigarOp.HARDMASK,
            "S": CigarOp.SOFTMASK
        }
        return

    def normalize(self, op_code):
        return self.known_ops[op_code]


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--paf-file", "-in",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="paf_file",
        help="PAF contig-to-contig alignment file (normalized)."
    )

    parser.add_argument(
        "--genome-size",
        "--fasta-index",
        "-s", "-f", "-gs",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="genome_size",
        default=None,
        help="Genome size file or FASTA index of the alignment target genome."
    )

    parser.add_argument(
        "--precision", "-prec",
        type=str,
        choices=["exact", "struct"],
        default="exact",
        dest="precision",
        help="Set desired precision: exact OR struct[ural]"
    )

    parser.add_argument(
        "--mapq-lower-bound", "-mq",
        type=int,
        default=0,
        dest="mapq_lower_bound",
        help=(
            "Retain only alignments that have a MAPQ above (>) "
            "this lower bound. Default: 0 "
            "(set to -1 to deactivate)"
        )
    )

    parser.add_argument(
        "--size-lower-bound", "-sz",
        type=int,
        default=0,
        dest="size_lower_bound",
        help=(
            "Retain only alignments where both the query and the target "
            "sequence have a length above (>) this lower bound. Default: 0"
        )
    )

    parser.add_argument(
        "--out-seq-stats", "-oss",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_seq_stats",
        default=None
    )

    parser.add_argument(
        "--verbose", "-vb",
        action="store_true",
        default=False,
        dest="verbose",
        help="Write process log to stdout. Default: False"
    )

    args = parser.parse_args()

    return args


def compute_qv(num_errors, ref_size):
    """_summary_

    Args:
        num_errors (_type_): _description_
        ref_size (_type_): _description_

    Returns:
        _type_: _description_
    """
    p = num_errors / ref_size
    try:
        q = -10 * math.log10(p)
    except ValueError:
        return 99
    return int(round(q, 0))


def read_genome_sizes(file_path, size_lower_bound):

    genome_sizes = dict()
    discard_log = []
    with open(file_path, "r") as listing:
        for line in listing:
            columns = line.strip().split()
            seq_name, seq_length = columns[:2]
            seq_length = int(seq_length)
            if seq_length <= size_lower_bound:
                discard_log.append(
                    (
                        "read_genome_sizes\tdiscard\t"
                        f"{seq_name}\t{seq_length}<={size_lower_bound}"
                    )
                )
            genome_sizes[seq_name] = int(seq_length)
    return genome_sizes, discard_log


def filter_alignments(alignments, mapq_bound, size_bound):

    discard_log = []
    discard = alignments["mapq"] <= mapq_bound
    if discard.sum() > 0:
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{discard.sum()}-aln-rows\tMAPQ<={mapq_bound}"
            )
        )
        # not sure what to record here ...
        #mapq_subset = alignments.loc[discard, :].copy()
        alignments = alignments.loc[~discard, :].copy()
    discard = (alignments["query_length"] <= size_bound) | (alignments["target_length"] <= size_bound)
    if discard.sum() > 0:
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{discard.sum()}-aln-rows\t[QRY|TRG]LEN<={size_bound}"
            )
        )
        size_subset = alignments.loc[discard, :].copy()
        queries_dropped = size_subset["query_name"].nunique()
        targets_dropped = size_subset["target_name"].nunique()
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{queries_dropped}-num-seq\tQRYLEN<={size_bound}"
            )
        )
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{targets_dropped}-num-seq\tTRGLEN<={size_bound}"
            )
        )
        alignments = alignments.loc[~discard, :].copy()

    return alignments, discard_log


def get_cigar_op_sets(precision):

    if precision == "exact":
        good_ops = set([CigarOp.IDENTITY])
        bad_ops = set([CigarOp.MISMATCH, CigarOp.DELETION])
        static_ops = set([CigarOp.INSERTION])
    elif precision == "struct":
        good_ops = set([CigarOp.IDENTITY, CigarOp.MISMATCH, CigarOp.MATCH])
        bad_ops = set([CigarOp.DELETION])
        static_ops = set([CigarOp.INSERTION])
    else:
        raise ValueError(f"Unknown level for precision parameter: {precision}")
    return good_ops, bad_ops, static_ops


def split_cigar_string(normalizer, size_dist, cigar_re, cigar_ops):

    op_series = []
    for cigar_op in cigar_re.finditer(cigar_ops):
        step = int(cigar_op.group("NUM"))
        norm_op = normalizer.normalize(cigar_op.group("OP"))
        op_series.append((step, norm_op))
        size_dist[norm_op].append(step)
    return op_series


def walk_cigar_ops(good_ops, bad_ops, static_ops, aln_start, aln_end, op_series):

    start = aln_start
    end = aln_start
    last_op = None
    for step, op_type in op_series:
        if op_type in good_ops:
            end += step
            last_op = op_type
        elif op_type in bad_ops:
            yield start, end
            # this steps over the
            # unsupported region
            start = end + step
            end = start
            last_op = op_type
        elif op_type in static_ops:
            last_op = op_type
        else:
            raise ValueError(f"Unsupported CIGAR operation: {op_type}")
    assert last_op in good_ops
    assert end == aln_end
    yield start, end

    return


#def compute_block_qv()


def process_alignment_file(alignments, precision):

    CigNorm = CigarNormalizer()
    cigar_op_sizes = col.defaultdict(list)
    CIGAR_OP_RE = re.compile("(?P<NUM>[0-9]+)(?P<OP>[MX\=DIHS]{1})")

    supported_regions = col.defaultdict(list)

    good_ops, bad_ops, static_ops = get_cigar_op_sets(precision)

    split_cigar = fnt.partial(split_cigar_string, CigNorm, cigar_op_sizes, CIGAR_OP_RE)
    walk_cigar = fnt.partial(walk_cigar_ops, good_ops, bad_ops, static_ops)

    for row in alignments.itertuples():
        op_series = split_cigar(row.cg_cigar)
        regions = [(start, end) for start, end in walk_cigar(row.target_start, row.target_end, op_series)]
        supported_regions[row.target_name].extend(regions)

    return cigar_op_sizes, supported_regions


def compute_sequence_statistics(seq_sizes, supported_regions):


    seq_stats = []
    for seq_name, seq_size in seq_sizes.items():
        try:
            seq_regions = supported_regions[seq_name]
            num_regions = len(seq_regions)
            total_support = sum(reg[1]-reg[0] for reg in seq_regions)
            total_error = seq_size - total_support
            qv_est = compute_qv(total_error, seq_size)
            support_pct = round(total_support/seq_size * 100, 2)
            seq_stats.append(
                (seq_name, seq_size, num_regions, total_support, support_pct, qv_est)
            )
        except KeyError:
            seq_stats.append((seq_name, seq_size, 0, 0, 0, 0))

    seq_stats = pd.DataFrame.from_records(
        seq_stats, columns=[
            "seq_name", "seq_size", "support_regions",
            "support_length", "support_pct", "support_qv"
        ]
    )

    seq_stats.sort_values(["seq_name"], inplace=True)

    return seq_stats


def log_process(log_lines, verbose):

    if log_lines and verbose:
        msg = "\n".join(log_lines) + "\n"
        sys.stdout.write(msg)
    return


def main():

    args = parse_command_line()

    gsize, log_lines = read_genome_sizes(args.genome_size)
    log_process(log_lines, args.verbose)

    paf_aln = pd.read_csv(args.paf_file, sep="\t", header=0)
    paf_aln = paf_aln.loc[paf_aln["tp_align_type"] != 2, :].copy()

    paf_aln, log_lines = filter_alignments(paf_aln, args.mapq_lower_bound, args.size_lower_bound)
    log_process(log_lines, args.verbose)

    _, supported_regions = process_alignment_file(paf_aln, args.precision)

    seq_stats = compute_sequence_statistics(gsize, supported_regions)

    if args.out_seq_stats is not None:
        args.out_seq_stats.parent.mkdir(parents=True, exist_ok=True)
        seq_stats.to_csv(args.out_seq_stats, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()

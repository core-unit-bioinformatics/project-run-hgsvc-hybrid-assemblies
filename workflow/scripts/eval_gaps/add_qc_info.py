#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl
import functools

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--gap-file", "-g",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="gap_file"
    )

    parser.add_argument(
        "--block-intersect", "-b",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="block_intersect"
    )

    parser.add_argument(
        "--nucfreq-regions", "-n",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="nucfreq_regions"
    )

    parser.add_argument(
        "--flagger-regions", "-f",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="flagger_regions"
    )

    parser.add_argument(
        "--chrom-assign", "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="chrom_assign"
    )

    args = parser.parse_args()

    return args


def load_gap_ids(file_path):

    with open(file_path) as listing:
        header = listing.readline().strip().split()
        assert header[0].startswith("#")
        header[0] = "chrom"
        assert header[3] == "name"
        header[3] = "gap_id"
        header[1] = "gap_start"
        header[2] = "gap_end"

    gaps = pd.read_csv(
        file_path, comment="#", sep="\t",
        header=None, names=header,
        usecols=header[:-2]
    )
    gaps = gaps[["chrom", "gap_start", "gap_end", "gap_id"]].copy()

    count_stats = {
        "total_gaps": gaps["gap_id"].nunique(),
        "autosome_gaps": (~gaps["chrom"].isin(["chrX", "chrY"])).sum(),
        "female_gaps": (gaps["chrom"].isin(["chrX"])).sum(),
        "male_gaps": (gaps["chrom"].isin(["chrY"])).sum(),
    }

    return gaps, count_stats


def load_chrom_assign(file_path):

    df = pd.read_csv(file_path, sep="\t", comment="#", header=0)

    chrom_assign = dict()

    for query, assignments in df.groupby("query_name"):
        max_match = assignments["align_matching"].idxmax()
        chrom = assignments.at[max_match, "target_name"]
        chrom_assign[query] = chrom

    return chrom_assign


def read_flagger_regions(file_path):

    df = pd.read_csv(file_path, sep="\t", header=0)
    df.rename({"#contig": "contig"}, axis=1, inplace=True)
    df.drop(["score", "asm_unit"], axis=1, inplace=True)
    return df


def read_nucfreq_regions(file_path):

    df = pd.read_csv(file_path, sep="\t", header=0)
    df.rename({"#contig": "contig"}, axis=1, inplace=True)
    df.drop(["hifi_median_cov", "hifi_pct_median_cov", "ont_median_cov", "ont_pct_median_cov"], axis=1, inplace=True)
    return df


def summarize_regions(regions, contig, start, end):

    select_ctg = regions["contig"] == contig
    select_start = regions["start"] < end
    select_end = regions["end"] > start

    selector = select_ctg & select_start & select_end
    if not selector.any():
        summary = [0, 0]
    else:
        sub = regions.loc[selector, :].copy()
        first = sub.index[0]
        last = sub.index[-1]
        sub.loc[first, "start"] = max(start, sub.at[first, "start"])
        assert sub.at[first, "end"] > start
        sub.loc[last, "end"] = min(end, sub.at[last, "end"])
        assert sub.at[last, "start"] < end
        total = (sub["end"] - sub["start"]).sum()
        pct = round(total / (end - start) * 100, 5)
        summary = [total, pct]
    return summary


@functools.lru_cache(2048)
def split_aln_context(align_context):

    _, aln_context = align_context.split("|")
    asm_ctg, asm_coord = aln_context.split(":")
    asm_start, asm_end = asm_coord.split("-")
    asm_start = int(asm_start)
    asm_end = int(asm_end)
    return asm_ctg, asm_start, asm_end



def process_gap_alnblock_table(file_path, flagger_file, nucfreq_file, gaps, chrom_assign):

    df = pd.read_csv(file_path, sep="\t", header=0)
    file_parts = file_path.name.split(".")
    sample = file_parts[0] + "." + file_parts[1]
    asm_unit = file_parts[2]
    df.drop(["aln_base_block", "aln_coarse_block"], axis=1, inplace=True)

    flagger = read_flagger_regions(flagger_file)
    nucfreq = read_nucfreq_regions(nucfreq_file)

    gap_status = []

    for row in gaps.itertuples():
        gap_id = row.gap_id
        gap_chrom = row.chrom
        sub = df.loc[df["gap_id"] == gap_id, :]
        if sub.empty or (sub["aln_label"] == "ASMGAP").all():
            status_info = {
                "gap_id": gap_id,
                "sample": sample, "asm_unit": asm_unit,
                "seq": "unknown", "start": -1, "end": -1,
                "overlap_pct": 0,
                "align_status": "open",
                "flagger_bp": 0, "flagger_pct": 0,
                "nucfreq_bp": 0, "nucfreq_pct": 0
            }
            gap_status.append(status_info)
        elif (sub["aln_label"] == "ALN").all() and (sub["overlap_pct"] > 99.9).all():

            match_found = False
            for row in sub.itertuples():
                ctg, start, end = split_aln_context(row.aln_info)
                ctg_chrom = chrom_assign[ctg]
                if ctg_chrom != gap_chrom:
                    continue
                flagger_summary = summarize_regions(flagger, ctg, start, end)
                nucfreq_summary = summarize_regions(nucfreq, ctg, start, end)
                status_info = {
                    "gap_id": gap_id,
                    "sample": sample, "asm_unit": asm_unit,
                    "seq": ctg, "start": start, "end": end,
                    "overlap_pct": row.overlap_pct,
                    "align_status": "covered",
                    "flagger_bp": flagger_summary[0], "flagger_pct": flagger_summary[1],
                    "nucfreq_bp": nucfreq_summary[0], "nucfreq_pct": nucfreq_summary[1]
                }
                gap_status.append(status_info)
                match_found = True

            if not match_found:
                status_info = {
                    "gap_id": gap_id,
                    "sample": sample, "asm_unit": asm_unit,
                    "seq": "unknown", "start": -1, "end": -1,
                    "overlap_pct": 0,
                    "align_status": "open",
                    "flagger_bp": 0, "flagger_pct": 0,
                    "nucfreq_bp": 0, "nucfreq_pct": 0
                }
                gap_status.append(status_info)

        else:
            if (sub["aln_label"] == "ALN").any():
                match_found = False
                for row in sub.itertuples():
                    if row.aln_label != "ALN":
                        continue
                    ctg, start, end = split_aln_context(row.aln_info)
                    ctg_chrom = chrom_assign[ctg]
                    if ctg_chrom != gap_chrom:
                        continue
                    flagger_summary = summarize_regions(flagger, ctg, start, end)
                    nucfreq_summary = summarize_regions(nucfreq, ctg, start, end)
                    status_label = "partial"
                    if row.overlap_pct > 99.9:
                        status_label = "covered"
                    status_info = {
                        "gap_id": gap_id,
                        "sample": sample, "asm_unit": asm_unit,
                        "seq": ctg, "start": start, "end": end,
                        "overlap_pct": row.overlap_pct,
                        "align_status": status_label,
                        "flagger_bp": flagger_summary[0], "flagger_pct": flagger_summary[1],
                        "nucfreq_bp": nucfreq_summary[0], "nucfreq_pct": nucfreq_summary[1]
                    }
                    gap_status.append(status_info)
                    match_found = True

                if not match_found:
                    status_info = {
                        "gap_id": gap_id,
                        "sample": sample, "asm_unit": asm_unit,
                        "seq": "unknown", "start": -1, "end": -1,
                        "overlap_pct": 0,
                        "align_status": "open",
                        "flagger_bp": 0, "flagger_pct": 0,
                        "nucfreq_bp": 0, "nucfreq_pct": 0
                    }
                    gap_status.append(status_info)

            else:
                status_info = {
                    "gap_id": gap_id,
                    "sample": sample, "asm_unit": asm_unit,
                    "seq": "unknown", "start": -1, "end": -1,
                    "overlap_pct": sub["overlap_pct"].sum(),
                    "align_status": "complex",
                    "flagger_bp": 0, "flagger_pct": 0,
                    "nucfreq_bp": 0, "nucfreq_pct": 0
                }
                gap_status.append(status_info)

    gap_status = pd.DataFrame.from_records(gap_status)
    return gap_status


def summarize_status(gap_table):

    thresholds = [101, 10, 1, 0.1]
    labels = ["closed_at_any", "closed_at_pct-10", "closed_at_pct-1", "closed_at_pct-01"]
    assert len(thresholds) == len(labels)

    summary = []
    stats = col.Counter()

    for gap_id, gap_status in gap_table.groupby("gap_id"):
        if gap_status.shape[0] > 1:
            if gap_status["align_status"].nunique() > 1:
                raise RuntimeError(gap_status)
            seqs = ";".join(sorted(gap_status["seq"].values))
            aln_label = gap_status["align_status"].iloc[0]
            if aln_label == "covered":
                raise RuntimeError(gap_status)
            summary.append(
                (
                    gap_id, gap_status["sample"].iloc[0], gap_status["asm_unit"].iloc[0],
                    seqs, aln_label, "not_closed"
                )
            )
            stats[aln_label] += 1
        else:
            if gap_status["align_status"].iloc[0] == "covered":
                flagger_pct = gap_status["flagger_pct"].iloc[0]
                nucfreq_pct = gap_status["nucfreq_pct"].iloc[0]

                last_label = None
                for t, l in zip(thresholds, labels):
                    if flagger_pct < t and nucfreq_pct < t:
                        stats[l] += 1
                        last_label = l

                summary.append(
                    (
                        gap_id, gap_status["sample"].iloc[0], gap_status["asm_unit"].iloc[0],
                        gap_status["seq"].iloc[0], gap_status["align_status"].iloc[0], last_label
                    )
                )

            else:
                summary.append(
                (
                    gap_id, gap_status["sample"].iloc[0], gap_status["asm_unit"].iloc[0],
                    gap_status["seq"].iloc[0], gap_status["align_status"].iloc[0], "not_closed"
                )
            )

    summary = pd.DataFrame.from_records(
        summary, columns=["gap_id", "sample", "asm_unit", "asm_seq", "aln_status", "stringency"]
    )

    return summary


def main():

    args = parse_command_line()

    gaps, gap_stats = load_gap_ids(args.gap_file)
    chrom_assign = load_chrom_assign(args.chrom_assign)

    gap_table = process_gap_alnblock_table(
        args.block_intersect,
        args.flagger_regions,
        args.nucfreq_regions,
        gaps,
        chrom_assign
    )

    summary = summarize_status(gap_table)
    summary = gaps.merge(summary, on="gap_id", how="outer")
    count_stats = summary.groupby(["chrom", "aln_status", "stringency"]).size().reset_index()
    count_stats.rename({0: "count"}, inplace=True)
    print(count_stats)
    # count_stats.update(gap_stats)
    # print(count_stats)
    # print(summary.head())

    gaps = gaps.merge(gap_table, on="gap_id")



    return 0


if __name__ == "__main__":
    main()

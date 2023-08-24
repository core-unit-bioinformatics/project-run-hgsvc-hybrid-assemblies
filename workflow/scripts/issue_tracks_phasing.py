#!/usr/bin/env python3

import argparse as argp
import collections as col
import hashlib as hl
import pathlib as pl

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--contig-cov",
        "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="contig_cov",
    )

    parser.add_argument(
        "--aln-unassign",
        "-u",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="aln_unassign",
    )

    parser.add_argument(
        "--aln-hap1",
        "-h1",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="aln_hap1",
    )

    parser.add_argument(
        "--aln-hap2",
        "-h2",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="aln_hap2",
    )

    parser.add_argument(
        "--acrocentrics",
        "-a",
        type=str,
        nargs="+",
        dest="acrocentrics",
        default=["chr13", "chr14", "chr15", "chr21", "chr22"]
    )

    args = parser.parse_args()

    return args


def load_contig_coverages(norm_paf, acrocentrics):

    cov = pd.read_csv(
        norm_paf, sep="\t",
        header=[0, 1, 2],
        index_col=[0, 1, 2]
    )

    # drop chrM, chrX, chrY (for now)
    select_not_auto = cov.index.get_level_values("chrom").isin(["chrM", "chrX", "chrY"])
    cov = cov.loc[~select_not_auto, :].copy()

    acro_rows = cov.index.get_level_values("chrom").isin(acrocentrics)

    rdna_has_cov = (cov[
        [("asm-rdna", "0", "ctg_align_cov"), ("asm-rdna", "60", "ctg_align_cov")]
    ] > 0).any(axis=1)

    # create a boolean mask that indicates rows in the
    # coverage table that are likely uncovered by
    # hap1/hap2/unassign because they represent rDNA
    # arrays (separate in Verkko)
    rdna_mask = rdna_has_cov & acro_rows
    rdna_none = ~rdna_mask

    return cov, rdna_none


def load_alignments(hap1, hap2, unassign):

    drop_columns = [
        "cg_cigar", "cs_tag_diffs", "s1_chaining_score", "s2_chaining_score",
        "rl_len_rep_seed_query", "cm_num_chain_miniz", "ms_dp_align_segment_max_score"
    ]

    aln_hap1 = pd.read_csv(hap1, sep="\t", header=0)
    aln_hap1.drop(drop_columns, axis=1, inplace=True)
    aln_hap1 = aln_hap1.loc[aln_hap1["tp_align_type"] == 1, :].copy()

    aln_hap2 = pd.read_csv(hap2, sep="\t", header=0)
    aln_hap2.drop(drop_columns, axis=1, inplace=True)
    aln_hap2 = aln_hap2.loc[aln_hap2["tp_align_type"] == 1, :].copy()

    aln_unassign = pd.read_csv(unassign, sep="\t", header=0)
    aln_unassign.drop(drop_columns, axis=1, inplace=True)
    aln_unassign = aln_unassign.loc[aln_unassign["tp_align_type"] == 1, :].copy()

    alignments = {
        "hap1": aln_hap1,
        "hap2": aln_hap2,
        "unassigned": aln_unassign
    }

    return alignments


def get_cov_selector(asm_unit):
    selector = [
        (f"asm-{asm_unit}", "0", "ctg_align_cov"),
        (f"asm-{asm_unit}", "60", "ctg_align_cov"),
    ]
    return selector


def create_target_regions(tview_regions, tview_aln_records, issue_label, min_size=0):

    plain_regions = []
    for block_name, region_list in tview_regions.items():
        chrom = region_list[0][0]
        start = min(t[1] for t in region_list)
        end = max(t[2] for t in region_list)
        assert start < end
        if (end - start) < min_size:
            continue
        aln_records = tview_aln_records[block_name]
        plain_regions.append(
            (chrom, start, end, block_name, issue_label, aln_records)
        )
    plain_regions = pd.DataFrame.from_records(
        plain_regions, columns=["chrom", "start", "end", "name", "label", "align_records"]
    )
    plain_regions.sort_values(["chrom", "start", "end"], inplace=True)

    return plain_regions


def create_query_regions(qview_nested_regions, qview_aln_records, issue_label, match_blocks):

    plain_regions = []
    for block_name, region_lists in qview_nested_regions.items():
        if block_name not in match_blocks:
            continue
        for query_name, query_coords in region_lists.items():
            start = min(t[0] for t in query_coords)
            end = max(t[1] for t in query_coords)
            assert start < end
            aln_records = qview_aln_records[block_name]
            plain_regions.append(
                (query_name, start, end, block_name, issue_label, aln_records)
            )
    plain_regions = pd.DataFrame.from_records(
        plain_regions, columns=["contig", "start", "end", "name", "label", "align_records"]
    )
    plain_regions.sort_values(["contig", "start", "end"], inplace=True)

    return plain_regions


def extract_diploid_regions(ctg_cov, rdna_free, alignments):

    select_h1 = get_cov_selector("hap1")
    select_h2 = get_cov_selector("hap2")
    select_un = get_cov_selector("unassigned")

    no_unassign = (ctg_cov[select_un] == 0).all(axis=1)
    any_h1 = (ctg_cov[select_h1] == 1).any(axis=1)
    any_h2 = (ctg_cov[select_h2] == 1).any(axis=1)

    selector = rdna_free & no_unassign & any_h1 & any_h2

    diploid_regions = ctg_cov.loc[selector, :]
    use_alignments = ["hap1", "hap2"]
    label = "diploid"

    tview_regions, qview_regions = process_selected_regions(
        diploid_regions, use_alignments, alignments, label
    )

    return tview_regions, qview_regions


def extract_dip_phasing_error(ctg_cov, rdna_free, alignments):

    select_h1 = get_cov_selector("hap1")
    select_h2 = get_cov_selector("hap2")
    select_un = get_cov_selector("unassigned")

    unassign = (ctg_cov[select_un] > 1).any(axis=1)
    no_h1 = (ctg_cov[select_h1] == 0).all(axis=1)
    no_h2 = (ctg_cov[select_h2] == 0).all(axis=1)

    selector = rdna_free & unassign & no_h1 & no_h2
    error_regions = ctg_cov.loc[selector, :]
    label = "dip_phasing_error"
    use_alignments = ["unassigned"]

    tview_regions, qview_regions = process_selected_regions(
        error_regions, use_alignments, alignments, label
    )

    return tview_regions, qview_regions


def extract_hap_phasing_error(ctg_cov, rdna_free, alignments, main, other):

    select_main = get_cov_selector(main)
    select_other = get_cov_selector(other)
    select_un = get_cov_selector("unassigned")

    unassign = (ctg_cov[select_un] > 0).any(axis=1)
    no_main = (ctg_cov[select_main] == 0).all(axis=1)
    some_other = (ctg_cov[select_other] > 0).any(axis=1)

    selector = rdna_free & unassign & no_main & some_other
    error_regions = ctg_cov.loc[selector, :]
    label = f"{main}_phasing_error"
    use_alignments = ["unassigned"]

    tview_regions, qview_regions = process_selected_regions(
        error_regions, use_alignments, alignments, label
    )

    return tview_regions, qview_regions


def extract_loh_regions(ctg_cov, rdna_free, alignments):

    select_h1 = get_cov_selector("hap1")
    select_h2 = get_cov_selector("hap2")
    select_un = get_cov_selector("unassigned")

    no_unassign = (ctg_cov[select_un] == 0).all(axis=1)
    hap_h1 = (ctg_cov[select_h1] == 1).all(axis=1)
    hap_h2 = (ctg_cov[select_h2] == 1).all(axis=1)

    selector = no_unassign & rdna_free & (hap_h1 ^ hap_h2)
    error_regions = ctg_cov.loc[selector, :]
    label = "LOH"
    use_alignments = ["hap1", "hap2"]

    tview_regions, qview_regions = process_selected_regions(
        error_regions, use_alignments, alignments, label, int(1e6)
    )

    return tview_regions, qview_regions


def process_selected_regions(ctg_cov, use_alignments, alignments, label, min_size=0):

    annotated_regions = annotate_regions_with_alignments(ctg_cov, use_alignments, alignments)
    tview_regions = create_target_regions(
        annotated_regions["target_view_regions"],
        annotated_regions["target_view_aln_records"],
        label, min_size
    )

    qview_regions = create_query_regions(
        annotated_regions["query_view_regions"],
        annotated_regions["query_view_aln_records"],
        label, set(tview_regions["name"].unique())
    )

    return tview_regions, qview_regions


def annotate_regions_with_alignments(issue_regions, use_alignments, alignments):

    # target view = BED-like track in reference coordinates
    target_view_regions = col.defaultdict(list)
    target_view_aln_records = dict()

    # query view = BED-like track in assembly coordinates
    query_view_regions = dict()
    query_view_aln_records = dict()

    for chrom, start, end in issue_regions.index:

        target_view_records = set()
        query_view_records = set()
        regions_by_query = col.defaultdict(list)
        for aln_name in use_alignments:
            aln = alignments[aln_name]
            select_chrom = aln["target_name"] == chrom
            select_start = aln["target_end"] >= start
            select_end = aln["target_start"] < end

            sub = aln.loc[select_chrom & select_start & select_end, :]
            # TODO: the below can be made more succinct
            # by preprocessing the alignments
            for row in sub.itertuples():
                # target/ref perspective
                # - what query aligns here?
                target_view_aln = (
                    f"QN:{row.query_name},"
                    f"QL:{row.query_length},"
                    f"OR:{row.align_orient},"
                    f"QS:{row.query_start},"
                    f"QE:{row.query_end},"
                    f"MQ:{row.mapq}"
                )
                target_view_records.add(target_view_aln)

                # query/assembly perspective
                # - where do I align to?
                query_view_aln = (
                    f"TN:{row.target_name},"
                    f"QL:{row.query_length},"
                    f"OR:{row.align_orient},"
                    f"TS:{row.target_start},"
                    f"TE:{row.target_end},"
                    f"MQ:{row.mapq}"
                )
                query_view_records.add(query_view_aln)
                regions_by_query[row.query_name].append((row.query_start, row.query_end))

        target_view_records = "|".join(sorted(target_view_records))
        query_view_records = "|".join(sorted(query_view_records))
        joined_name = target_view_records + "@" + query_view_records
        block_name = hl.md5(joined_name.encode("utf-8")).hexdigest()

        target_view_regions[block_name].append((chrom, start, end))
        target_view_aln_records[block_name] = target_view_records

        query_view_regions[block_name] = regions_by_query
        query_view_aln_records[block_name] = query_view_records

    annotated = {
        "target_view_regions": target_view_regions,
        "target_view_aln_records": target_view_aln_records,
        "query_view_regions": query_view_regions,
        "query_view_aln_records": query_view_aln_records
    }
    return annotated



def main():

    args = parse_command_line()

    cov, rdna_free = load_contig_coverages(args.contig_cov, args.acrocentrics)

    alignments = load_alignments(args.aln_hap1, args.aln_hap2, args.aln_unassign)

    #dip_regions = extract_diploid_regions(cov, rdna_free, alignments)
    #dip_phase_error = extract_dip_phasing_error(cov, rdna_free, alignments)

    #hap1_phase_error = extract_hap_phasing_error(cov, rdna_free, alignments, "hap1", "hap2")
    #hap2_phase_error = extract_hap_phasing_error(cov, rdna_free, alignments, "hap2", "hap1")

    loh_regions = extract_loh_regions(cov, rdna_free, alignments)

    t, q = loh_regions

    print(t)
    print('===')
    print(q)

    return


if __name__ == "__main__":
    main()

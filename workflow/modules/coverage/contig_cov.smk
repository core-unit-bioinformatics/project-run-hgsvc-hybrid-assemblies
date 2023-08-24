
rule verkko_dump_assembly_coverage_issue_table:
    input:
        tsv = WORKDIR_EVAL.joinpath(
            "results/coverage/contig_ref/t2tv2",
            "{sample}.t2tv2.win-ctg-cov.tsv.gz"
        )
    output:
        bed = DIR_RES.joinpath(
            "issues", "contig_ref",
            "{sample}.t2tv2.ctg-cov-issues.bed.gz"
        )
    run:
        import pandas as pd

        select_cols = dict()
        for asm_unit in ["hap1", "hap2", "unassigned", "disconnected", "rdna"]:
            select_cov = [
                (f"asm-{asm_unit}", "MQ00", "ctg_align_cov"), (f"asm-{asm_unit}", "MQ60", "ctg_align_cov")
            ]
            select_cols[asm_unit] = select_cov
        select_cols["dip"] = select_cols["hap1"] + select_cols["hap2"]

        wg = pd.read_csv(input.tsv, sep="\t", header=[0, 1, 2], index_col=[0, 1, 2])
        # drop chrM
        not_mito = wg.index.get_level_values("chrom") != "chrM"
        wg = wg.loc[not_mito, :].copy()

        # build a boolean mask

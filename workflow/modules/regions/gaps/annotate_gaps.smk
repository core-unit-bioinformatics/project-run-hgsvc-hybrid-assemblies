
rule segdups_to_bed:
    input:
        sd95 = rules.split_segdup_annotation.output.sd95,
        sd98 = rules.split_segdup_annotation.output.sd98
    output:
        bed_95 = temp(DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-095.bed"
        )),
        bed_98 = temp(DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-098.bed"
        ))
    run:
        import pandas as pd

        for input_file, output_file, region_label in zip(input, output, ["SD95", "SD98"]):
            df = pd.read_csv(input_file, sep="\t", header=0)
            df.rename({"seq": "#contig"}, axis=1, inplace=True)
            df["name"] = region_label
            df.sort_values(["#contig", "start", "end"], inplace=True)
            df.to_csv(output_file, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_overlapping_segdups:
    input:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-{pct_id}.bed"
        ),
    output:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-{pct_id}.mrg.bed"
        ),
    wildcard_constraints:
        pct_id = "(095|098)"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    shell:
        "bedtools merge -c 4 -o distinct -i {input.bed} > {output.bed}"


rule annotate_gaps_with_segdups:
    """Annotate assemblies (assembly coordinate space)
    with segdup annotation. The input file in ref coordinates
    is just here for enforcing.
    """
    input:
        qry_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.asm-coord.bed"
        ),
        trg_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.ref-coord.bed"
        ),
        segdups = rules.merge_overlapping_segdups.output.bed
    output:
        isect = DIR_PROC.joinpath(
            "regions", "gaps", "intersections",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.isect-sd{pct_id}.tsv"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wao -a {input.qry_view} -b {input.segdups} > {output.isect}"


localrules: hprc_gaps_to_bed
rule hprc_gaps_to_bed:
    input:
        tsv = HPRC_COMMON_GAPS
    output:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "hprc_common_gaps.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(tsv, sep="\t", header=0)
        df.sort_values(["chrom", "start", "end"], inplace=True)
        df.rename({"chrom": "#chrom"}, axis=1, inplace=True)
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule annotate_gaps_with_hprc_gaps:
    """NB: HPRC gaps only exist in T2T space
    """
    input:
        qry_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.asm-coord.bed"
        ),
        trg_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.ref-coord.bed"
        ),
        gaps = rules.hprc_gaps_to_bed.output.bed
    output:
        isect = DIR_PROC.joinpath(
            "regions", "gaps", "intersections",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.isect-hprcgaps.tsv"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wao -a {input.trg_view} -b {input.gaps} > {output.isect}"


rule run_all_annotate_gaps:
    input:
        tables = expand(
            rules.annotate_gaps_with_segdups.output.isect,
            sample=SAMPLES,
            asm_unit=MAIN_ASSEMBLY_UNITS,
            refgenome=["hg38", "t2tv2"],
            pct_id=["095", "098"]
        ),
        gaps = expand(
            rules.annotate_gaps_with_hprc_gaps.output.isect,
            sample=SAMPLES,
            asm_unit=MAIN_ASSEMBLY_UNITS,
            refgenome=["t2tv2"]
        )

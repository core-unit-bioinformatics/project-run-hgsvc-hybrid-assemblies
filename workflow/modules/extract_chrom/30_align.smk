
assert SAMPLE_SEX is not None, "Need sample sex info for tig alignment"


rule create_extract_chrom_seq_idx:
    input:
        fasta = rules.fetch_tigs_from_sequence_files.output.fasta,
    output:
        fai = DIR_PROC.joinpath(
            "extract_chrom", "fasta_seqs",
            "{sample}.{select_chrom}.oriented.fasta.fai"
        )
    conda:
        DIR_ENVS.joinpath("samtools.yaml")
    shell:
        "samtools faidx {input.fasta}"


rule minimap_chrom_to_ref_align:
    input:
        fasta = rules.fetch_tigs_from_sequence_files.output.fasta,
        fai = rules.create_extract_chrom_seq_idx.output.fai,
        ref = select_ref_seq
    output:
        paf = DIR_RES.joinpath(
            "extract_chrom", "chrom_to_ref",
            "minimap", "{sample}.{ref}.{chrom}.seq.paf"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    shell:
        "minimap2 -x asm20 -L --eqx --MD --cs -c --secondary=no "
        "-t {threads} {input.ref} {input.fasta}"
        " | "
        "pigz > {output.paf}"


rule minimap_seqclass_to_chrom_align:
    input:
        fasta = rules.fetch_tigs_from_sequence_files.output.fasta,
        fai = rules.create_extract_chrom_seq_idx.output.fai,
        seqcls = select_seqclasses
    output:
        paf = DIR_RES.joinpath(
            "extract_chrom", "seqclass_to_chrom",
            "minimap", "{sample}.{ref}.{chrom}.cls.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    shell:
        "minimap2 -x asm20 -L --eqx --MD --cs -c --secondary=yes -N 10 -p 0.95 "
        "-t {threads} {input.fasta} {input.seqcls}"
        " | "
        "pigz > {output.paf}"


rule run_all_minimap_alignments:
    input:
        paf_ref = expand(
            rules.minimap_chrom_to_ref_align.output.paf,
            sample=MALE_SAMPLES,
            chrom=["chrY"],
            ref=["hg38", "t2t"]
        ),
        paf_cls = expand(
            rules.minimap_seqclass_to_chrom_align.output.paf,
            sample=MALE_SAMPLES,
            chrom=["chrY"],
            ref=["hg38", "t2t"]
        )
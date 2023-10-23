
rule run_ribotin_on_verkko:
    input:
        verkko_asm_check = WORKDIR_ASSEMBLY.joinpath(
            "proc/10-assemble/verkko", "{sample}.ps-none.ok"
        )
    output:
        checkfile = DIR_PROC.joinpath(
            "rdna", "ribotin", "{sample}.ps-none.ribotin.ok"
        )
    conda:
        DIR_ENVS.joinpath("ribotin.yaml")
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 32768,
        time_hrs=lambda wildcards, attempt: attempt * 11
    params:
        asm_folder = lambda wildcards, input: pathlib.Path(input.verkko_asm_check).with_suffix(".wd"),
        out_folder = lambda wildcards, output: pathlib.Path(output.checkfile).with_suffix(".wd"),
    shell:
        "ribotin-verkko --sample-name {wildcards.sample} -x human "
        "-i {params.asm_folder} -o {params.out_folder} "
        "--mbg `which MBG` --graphaligner `which GraphAligner`"
            " && "
        "touch {output.checkfile}"


rule run_all_ribotin:
    input:
        checkfiles = expand(
            rules.run_ribotin_on_verkko.output.checkfile,
            sample=SAMPLES
        )

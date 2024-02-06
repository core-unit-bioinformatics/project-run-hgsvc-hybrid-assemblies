import pandas

SAMPLE_TABLE = pandas.read_csv(
    config["samples"],
    sep="\t",
    comment="#",
    header=0
)

SAMPLES = sorted(SAMPLE_TABLE["sample"].unique())

if all("vrk-ps" in s for s in SAMPLES):

    ASSEMBLER = "verkko"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2", "unassigned"]

elif all("hsm-ps" in s for s in SAMPLES):

    ASSEMBLER = "hifiasm"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2"]

else:
    raise ValueError("Can process either Verkko or hifiasm assemblies, not both.")

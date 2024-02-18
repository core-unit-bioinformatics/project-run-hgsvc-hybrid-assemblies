import pandas

SAMPLE_TABLE = pandas.read_csv(
    config["samples"],
    sep="\t",
    comment="#",
    header=0
)

SAMPLES = sorted(SAMPLE_TABLE["sample"].unique())

SAMPLE_SEX = None
MALE_SAMPLES = None
FEMALE_SAMPLES = None
if "sex" in SAMPLE_TABLE:
    SAMPLE_SEX = dict()
    for row in SAMPLE_TABLE.itertuples():
        SAMPLE_SEX[row.sample] = row.sex
    MALE_SAMPLES = [s for s in SAMPLES if SAMPLE_SEX[s] in ["male", "m"]]
    FEMALE_SAMPLES = [s for s in SAMPLES if SAMPLE_SEX[s] in ["female", "f"]]

if all("vrk-ps" in s for s in SAMPLES):

    ASSEMBLER = "verkko"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2", "unassigned"]

elif all("CEPH" in s for s in SAMPLES):

    # special setting for pedigree samples

    ASSEMBLER = "verkko"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2", "unassigned"]

elif all("hsm-ps" in s for s in SAMPLES):

    ASSEMBLER = "hifiasm"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2"]

else:
    raise ValueError("Can process either Verkko or hifiasm assemblies, not both.")

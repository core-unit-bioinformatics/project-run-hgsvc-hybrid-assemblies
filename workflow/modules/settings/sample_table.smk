import pandas

SAMPLE_TABLE = pandas.read_csv(
    config["samples"],
    sep="\t",
    comment="#",
    header=0
)

SAMPLES = sorted(SAMPLE_TABLE["sample"].unique())

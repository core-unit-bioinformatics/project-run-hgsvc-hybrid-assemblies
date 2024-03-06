import pandas


def load_assembly_unit_karyotypes(file_path):

    df = pandas.read_csv(file_path, sep="\t", comment="#", header=0)

    lut = dict()
    for sample, au_infos in df.groupby("sample"):
        sex_hap1 = au_infos["asm-hap1"].iloc[0]
        sex_hap2 = au_infos["asm-hap2"].iloc[0]
        sex_label = f"{sex_hap1}|{sex_hap2}"
        lut[sample] = sex_label
    return lut

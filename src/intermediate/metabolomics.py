import numpy as np
import openpyxl as opx
import pandas as pd
import gaussianprocessderivatives as gp
import os

basedir = os.path.expanduser("~/huo2025/")

datadir = basedir + "data/fig3/"

# load OD data
wb = opx.load_workbook(datadir + "OD_Check_20211105_fixed.xlsx")
for i in [1, 2]:
    wb.remove(wb["Sheet" + str(i)])
conditions = ["0.1% Gal, 0.4% Pal"] * 3 + ["0.1% Gal"] * 3
time = [0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90]
strains = ["365_1", "365_2", "365_3"]
od1 = []
od2 = []
datetime = []
cond = []
expt_time = []
strain = []
for sheet, t in zip(wb, time[::-1]):
    datetime.extend([sheet["B40"].value] * 7)
    cond.extend(conditions + ["Null"])
    expt_time.extend([t] * 7)
    strain.extend(strains * 2 + ["Null"])
    for row in sheet.iter_rows(min_row=29, max_row=35, min_col=2, max_col=3):
        od_raw = [c.value for c in row]
        od1.append(od_raw[0])
        od2.append(od_raw[1])
# initialise
results_dict = {
    "OD_1": od1,
    "OD_2": od2,
    "Datetime": datetime,
    "Condition": cond,
    "Time": expt_time,
    "Strain": strain,
}
df = pd.DataFrame(results_dict)
# Split the two technical measurement
df = pd.melt(
    df,
    id_vars=["Time", "Datetime", "Condition", "Strain"],
    value_vars=["OD_1", "OD_2"],
    var_name="Meas_rep",
    value_name="OD",
)
meas_rep_map = lambda x: x.split("_")[-1]
df["Meas_rep"] = df["Meas_rep"].apply(meas_rep_map)

# correct media
gdf = df.groupby(by=["Condition", "Strain", "Time"])["OD"].agg(["mean", "std"])
gdf.columns = ["OD mean", "OD err"]
gdf.reset_index(inplace=True)


def correct_media(row, df):
    filt_null = (df["Strain"] == "Null") & (df["Time"] == row["Time"])
    od_null = df.loc[filt_null, "OD mean"].to_numpy()[0]
    od_null_err = df.loc[filt_null, "OD err"].to_numpy()[0]
    row["Corrected OD"] = row["OD mean"] - od_null
    if row["Strain"] != "Null":
        # Propagate error
        row["Corrected OD err"] = np.sqrt(row["OD err"] ** 2 + od_null_err**2)
    else:  # We should not correct Null by Null
        row["Corrected OD err"] = 0.0
    return row


gdf = gdf.apply(correct_media, args=(gdf,), axis=1)

# correct OD
x = []
y = []
od_corr_file = datadir + "dilution_data_xiao.tsv"
with open(od_corr_file) as ocf:
    for line in ocf:
        x.append(float(line.split("\t")[0]))
        y.append(float(line.split("\t")[1]))
ocdf = pd.DataFrame({"x": x, "y": y})
ocdf.sort_values(by="x", inplace=True)
x = ocdf["x"].to_list()
y = ocdf["y"].to_list()
bd = {0: (-4, 4), 1: (-4, 4), 2: (-2.5, 0)}
g = gp.maternGP(bd, x, y)
g.findhyperparameters()
g.results()
g.predict(x, derivs=1)
x = gdf["Corrected OD"].to_list()
g.predict(x)
gdf["CC OD"] = pd.Series(g.f)
# No error propagation yet from OD err to GP err
gdf["GP err"] = pd.Series(g.fvar)

# metabolomics data
met = pd.read_excel(
    datadir + "Galactose and Palatinose Analysis.xlsx",
    sheet_name=2,
    header=1,
    skiprows=[2, 3],
    index_col=0,
    usecols=list(range(6)),
)
time = [0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90]
t_to_t = lambda x: time[int(x[1:]) - 1] if isinstance(x, str) else x
met["Time"] = met["Time"].apply(t_to_t)
biorep = lambda x: str(x).split(".")[-1] if not np.isnan(x) else x
met["Biorep"] = met["Strain"].apply(biorep)
met.drop("Delft media", axis="index", inplace=True)


def normalise_by_init(row, df):
    filt = (df["Time"] == 0.0) & (df["Biorep"] == row["Biorep"])
    init_row = df.loc[filt]
    row["Normalised Galactose"] = (row["Galactose"] / init_row["Galactose"]).to_numpy()[
        0
    ]
    row["Normalised Palatinose"] = (
        row["Palatinose"] / init_row["Palatinose"]
    ).to_numpy()[0]
    return row


met = met.apply(normalise_by_init, args=(met,), axis=1)
met = pd.melt(
    met,
    id_vars=["Strain", "Time", "Biorep"],
    value_vars=["Normalised Galactose", "Normalised Palatinose"],
    var_name="Sugar",
    value_name="Area",
)

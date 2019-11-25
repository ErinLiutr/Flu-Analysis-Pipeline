import pandas as pd
import os

os.chdir("../../Desktop/lab/flu_transmission/data/Reanalysis_of_EMIT/transmissionTree")
print os.getcwd()

# read dna data
data = {}
for f in os.listdir("fine/distinct"):
    s_id = f.split("_")[0]
    b_id = f.split("_")[3].split(".")[0]
    data[s_id] = b_id

# read metadata
df = pd.read_csv("metadata.csv")
rows = df.loc[df['JCVI_bac_id'].isin(data.values())]
date = dict(zip(rows["JCVI_bac_id"], rows["date_collection"]))

new_date = {}
for key, val in date.items():
    new_key = str(int(key))
    new_date[new_key] = val


# map id with collection date
c_date = {}
for id in data.keys():
    c_date[int(id)] = new_date[data[id]]

c = pd.Series(c_date)
c.to_csv("fine.csv")

import pandas as pd

df = pd.read_csv("/home/ubuntu/storage/second_run/thesis/ThermoBase.csv")

tk = pd.read_csv(
    "/home/ubuntu/storage/second_run/thesis/taxonkit_out.txt",
    sep="\t",
    names=["names", "taxid"]
)

# remove fake taxids (when taxid == name)
tk.loc[tk["names"] == tk["taxid"], "taxid"] = None

df = df.merge(tk, on="names", how="left")

df.to_csv("/home/ubuntu/storage/second_run/thesis/Thermobase_taxid.csv", index=False)

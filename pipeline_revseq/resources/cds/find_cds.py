import pandas as pd
import os

files=[ file for file in os.listdir() if "jsonl" in file ]
cds=[]
for file in files:
        name=file.split("_annotation")[0]
        print(name)
        anno=pd.read_json(path_or_buf=file, lines=True)
        anno2=anno[anno['accession'].str.contains(name)]
        cdsanno=anno2["genes"].tolist()
        for entry in cdsanno[0]:
                line={"name":name, "type":"cds", "start":entry["cds"][0]["nucleotide"]["range"][0]["begin"], "end":entry["cds"][0]["nucleotide"]["range"][0]["end"]}
                cds.append(line)

cds2 = pd.DataFrame(cds)
cds2.to_csv("cds.tsv", sep="\t", index=False)

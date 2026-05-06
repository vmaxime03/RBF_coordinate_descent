import json
import os
import pandas as pd
import glob 
import numpy as np 
import sys 

paths = sorted(glob.glob(os.path.join(sys.argv[1], "*sdf.json")))

frames = []
for path in paths:
    with open(path) as f:
        data = json.load(f)
    frames.append(pd.DataFrame(data["points"]))

all_data = pd.concat(frames, ignore_index=True)

all_data["bnorm"] = np.sqrt(all_data["bx"]**2 + all_data["by"]**2);

cols = ["px", "py", "alpha", "bx", "by", "s", "bnorm"]

print("PEARSON")
print(all_data[cols].corr("pearson").round(3))

print("SPEARMAN")
print(all_data[cols].corr("spearman").round(3))

print("KENDALL")
print(all_data[cols].corr("kendall").round(3))

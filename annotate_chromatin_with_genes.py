#! /usr/bin/env python3

"""This script will take a .bed file containing open chromatin regions, and a
.txt file containing annotated chromatin regions, and combines their
information.

This script requires the Python library tqdm"""

import sys
from tqdm import tqdm
import re

if len(sys.argv) < 2:
    print("Please provide a path to the .bed file")
    sys.exit(1)

if len(sys.argv) < 3:
    print("Please provide a path to the .txt annotation file")
    sys.exit(1)

if len(sys.argv) < 4:
    print("Please provide a location to save to")
    sys.exit(1)

KEYS = ["Nearest PromoterID"]

# Open files and get data
with open(sys.argv[1]) as f:
    region_data = f.read().splitlines()
with open(sys.argv[2]) as f:
    annotation_data = f.read().splitlines()

# Make dicts of region data
labels = ["Chr", "Start", "End", "TF", "Val", "Strand"]
regions = [{
 labels[i]: int(val) if 1 <= i <= 2 else float(val) if i == 4 else val\
  for i, val in enumerate(line.split())
} for line in region_data]

# Make dicts of annotation data
labels = annotation_data[0].split("\t")
labels[0] = labels[0][:6]
annotations = [{
 labels[i]: val for i, val in enumerate(line.split("\t"))
} for line in annotation_data[1:]]
for annotation in annotations:
    for key in annotation:
        try:
            annotation[key] = int(annotation[key])
        except: pass

# Make chromosomes collection
chromosomes = {f"chr{i}": {} for i in range(1, 23)}
chromosomes["chrX"] = {}
chromosomes["chrY"] = {}

# Identify overlaps in each chromosome
for chr_id, chromosome in chromosomes.items():
    print(chr_id)
    chromosome["regions"] = [region for region in regions if region["Chr"] == chr_id]
    chromosome["regions"].sort(key=lambda r: r["Start"])
    chromosome["annotations"] = [annotation for annotation in annotations if annotation["Chr"] == chr_id]
    chromosome["annotations"].sort(key=lambda a: a["Start"])

    # Find overlap
    for annotation in tqdm(chromosome["annotations"]):
        for region in chromosome["regions"]:
            if (annotation["Start"] <= region["Start"] and annotation["End"] >= region["Start"]) or (annotation["Start"] <= region["End"] and annotation["End"] >= region["End"]):
                for key in KEYS:
                    region[key] = annotation[key]
                break
    
    # Annotate regions with no overlap
    for region in chromosome["regions"]:
        for key in KEYS:
            region[key] = region.get(key) or "N/A"

# Output
lines = []
lines.append("\t".join(list(chromosomes["chr1"]["regions"][0].keys())))
for chromosome in chromosomes.values():
    lines += ["\t".join([str(val) for val in region.values()]) for region in chromosome["regions"]]
   
with open(f"{sys.argv[3]}/annotated_regions.txt", "w") as f:
    f.write("\n".join(lines))
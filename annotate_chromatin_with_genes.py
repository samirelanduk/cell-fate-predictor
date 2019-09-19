#! /usr/bin/env python3

"""This script will take a .bed file containing open chromatin regions, and a
.txt file containing annotated chromatin regions, and combines their
information."""

import sys
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

# Open files and get data
with open(sys.argv[1]) as f:
    region_data = f.read().splitlines()
with open(sys.argv[2]) as f:
    annotation_data = f.read().splitlines()

# Make dicts of region data
labels = ["Chr", "Start", "End", "Peak", "Val"]
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

# Check all IDs are unique
if len(set([region["Peak"] for region in regions])) != len(regions):
    print("Region Peak IDs are not unique!")
    sys.exit(1)
if len(set([annotation["PeakID"] for annotation in annotations])) != len(annotations):
    print("Annotation Peak IDs are not unique!")
    sys.exit(1)

# Make IDs keys
regions = {region["Peak"]: region for region in regions}
annotations = {annotation["PeakID"]: annotation for annotation in annotations}

# Combine data
combinations = []
for region in regions:
    combination = {**regions[region]}
    annotation = annotations[region]
    for key in ["Nearest PromoterID"]:
        combination[key] = annotation[key]
    combinations.append(combination)

# Output
lines = []
lines.append("\t".join(list(combinations[0].keys())))
lines += ["\t".join([str(val) for val in combo.values()]) for combo in combinations]
   
with open(f"{sys.argv[3]}/annotated_regions.txt", "w") as f:
    f.write("\n".join(lines))

    
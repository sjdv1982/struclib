"""
Executes snakegraph.seamless
First, the input files are bound
Then, the graph is equilibrated
Finally, all virtual files starting with "outputs/" are written 

For comparison, run snakemake on the original data
 (with --use-singularity, don't forget to adapt the oscar-star rule)
 This will also write the results in outputs/ (i.e. create a backup first)
Except for the output of oscar-star, which is non-deterministic,
 all results should be the same
"""
import seamless
from seamless.highlevel import Context
import json, os

cache = seamless.RedisCache()

print("Load graph...")
graph = json.load(open("snakegraph.seamless"))
ctx = seamless.highlevel.load_graph(graph)

print("Bind files...")
inputs = (
    "params/cluster-cutoff",
    "params/selected-cluster",
    "receptor.pdb",
    "ligand.pdb",
    "receptor-bound.pdb",
    "ligand-bound.pdb",
    "normal-modes.dat",
    "docking-result.dat",
    "docking-result-pairwise-lrmsd.txt"
)
for file in inputs:
    print(file)
    data = open(file).read()
    ctx.filesystem[file] = data

print("Equilibrate...")
ctx.equilibrate()
ctx.translate(force=True) ## kludge, is it necessary?
ctx.equilibrate()
fs = ctx.filesystem.value.data.value

print("Writing outputs:")
os.system("mkdir -p outputs")
print(list(fs.keys()))
for file in sorted(list(fs.keys())):    
    if not file.startswith("outputs"):
        continue
    value = fs[file]
    if value is None:
        continue
    print(file)
    with open(file, "w") as f:
        f.write(str(value))

    

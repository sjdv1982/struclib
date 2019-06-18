import requests, os

result = {}
for pdbcode in pdbcodes:
    outfile = "%s.pdb" % pdbcode
    if not os.path.exists(outfile) and len(pdbcode) == 4:
        pdb = requests.get("https://files.rcsb.org/download/%s.pdb" % pdbcode).text
        with open(outfile, "w") as f:
            f.write(pdb)
    else:
        pdb = open(outfile).read()
    result[pdbcode] = pdb
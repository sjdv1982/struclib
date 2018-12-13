from seamless.silk import Silk
import requests
StructureState = Silk(schema=structurestate_schema.data, data=None)

struc = StructureState()

for pdbcode in pdbcodes:
    pdb = requests.get("https://files.rcsb.org/download/%s.pdb" % pdbcode).text
    struc.parse_pdb(pdb, obj=pdbcode)
result = struc.data

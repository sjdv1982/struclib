from seamless.silk import Silk
import requests
StructureState = Silk(schema=structurestate_schema, data=None)

struc = StructureState()

for pdbcode, pdb in pdbs.items():
    struc.parse_pdb(pdb, obj=pdbcode)
result = struc.data
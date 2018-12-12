import os; os.chdir(os.path.dirname(__file__))
from struclib.StructureState import StructureState
pdb = open("1AVX.pdb").read()
s = StructureState()
s.parse_pdb(pdb)
s.select("resname == 'ALA' and model == 1 and resid < 200")
s.get_selection()

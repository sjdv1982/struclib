import os; os.chdir(os.path.abspath(os.path.dirname(__file__)))
from struclib.StructureState import StructureState
pdb = open("1AVX.pdb").read()
s0 = StructureState()
from seamless.silk import Silk
s = Silk(schema=s0.schema)()
s.parse_pdb(pdb, obj="1AVX")
s.select("resname == 'ALA' and model == 1 and resid < 200")
s.get_selection()
s.parse_pdb(pdb, obj="1AKE")

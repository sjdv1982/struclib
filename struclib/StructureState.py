from seamless.silk import Silk
from seamless.silk.meta import meta, validator

class StructureState(metaclass=meta):
    
    # DTypes of atomic data and atomic state
    atomic_dtype = [
        ("model", 'uint16'),            
        ("hetero", "S1"),
        ("name", "S4"),
        ("altloc","S1"),
        ("resname", "S3"),            
        ("chain","S1"),
        ("index", 'uint32'),
        ("icode", "S1"), 
        ("resid", 'uint16'),            
        ("x", 'float32'),
        ("y", 'float32'),
        ("z", 'float32'),
        ("occupancy", 'float32'),
        ("bfactor", 'float32'),
        ("segid", "S4"),
        ("element", "S2")                  
    ]
    atomic_state_dtype = [
        ("obj", 'uint8')] + atomic_dtype + \
        [("sele", 'uint8'), ("repr", 'uint16')
    ]

    def __init__(self):
        import numpy as np        
        atomic_state_dtype = np.dtype(self.atomic_state_dtype, align=True)

        # Array of the atomic state
        self.atomstate = np.zeros(0, dtype=atomic_state_dtype)
        
        # Names of molecular objects
        self.PRIVATE_objects = []  
        
        # Names of representations
        self.PRIVATE_repr = ["cartoon", "spacefill", "lines", "sticks", "ball+sticks", "dots"]  
        
        # Names of selections
        self.PRIVATE_sele = []
    
    def parse_pdb(self, pdb, obj=None):
        import Bio.PDB
        from io import StringIO
        import numpy as np
        
        # Parse from file object, or from string
        if callable(getattr(pdb, "read", None)):
            pdb_obj = pdb
        else:
            pdb_obj = StringIO(pdb)
        if obj is None:
            for count in range(256):
                obj = "obj%0.2d" % (count+1)
                if obj not in self.PRIVATE_objects:
                    break
        if obj in self.PRIVATE_objects:
            raise ValueError("object %s already exists" % obj)
        self.PRIVATE_objects.append(obj)
        
        p = Bio.PDB.PDBParser()
        struc = p.get_structure(obj, pdb_obj)
        natoms = len(list(struc.get_atoms()))
        new_statelen = len(self.atomstate) + natoms
        atomstate = np.zeros(new_statelen,dtype=self.atomstate.dtype)
        assert atomstate.dtype == self.atomstate.data.dtype
        atomstate[:len(self.atomstate)] = self.atomstate.data
        
        a = atomstate[len(self.atomstate):]
        a["repr"] = 1  # cartoon
        a["obj"] = len(self.PRIVATE_objects)
        count = 0
        for modelnr, model in enumerate(struc.get_models()):
            for chain in model.get_chains():
                chainid = chain.id
                for residue in chain.get_residues():  
                    hetero, resid, icode = residue.get_id()
                    segid = residue.segid
                    resname = residue.resname
                    for atom in residue.get_atoms():
                        aa = a[count]
                        aa["model"] = modelnr + 1
                        aa["hetero"] = hetero
                        aa["name"] = atom.name
                        aa["altloc"] = atom.altloc
                        aa["resname"] = resname
                        aa["chain"] = chainid
                        aa["index"] = atom.serial_number
                        aa["icode"] = icode
                        aa["resid"] = resid
                        aa["x"] = atom.coord[0]
                        aa["y"] = atom.coord[1]
                        aa["z"] = atom.coord[2]
                        occ = atom.occupancy
                        if occ is None or occ < 0:
                            occ = 0
                        aa["occupancy"] = occ
                        aa["segid"] = segid
                        aa["element"] = atom.element
                        count += 1
        self.atomstate = atomstate
        
    def select(self, query, sele="sele"):        
        import pandas as pd # Must be patched version of pandas!! spin off pandas.core.computation? "pandeval"?
        import numpy as np
        old_selearray = None
        sele_state = self.atomstate.data["sele"]
        try:
            selenr = self.PRIVATE_sele.index(sele) + 1
            old_selearray = (sele_state == selenr)  
            sele_state[old_selearray] = 0
        except ValueError:
            self.PRIVATE_sele.append(sele)
            selenr = len(self.PRIVATE_sele)
        dic = {key: self.atomstate.data[key] for key in self.atomstate.dtype.fields.keys() \
               if key not in ("obj", "repr")}
        objects = np.array([""] + self.PRIVATE_objects.data).astype("S")
        dic["all"] = np.ones(len(self.atomstate))
        dic["obj"] = objects[self.atomstate.data["obj"]]
        for xselenr, xsele in enumerate(self.PRIVATE_sele):
            if xsele in dic or (xselenr == selenr and old_selearray is None):
                continue
            xselearray = (sele_state == (xselenr + 1) )
            print("XSELE", xsele, xselearray.sum())
            dic[xsele] = xselearray
        try:
            selearray = pd.eval(query, global_dict=dic, align_result=False, str_as_bytes=True)
            print("%d atoms selected" % selearray.sum())
            sele_state[selearray] = selenr
        except Exception:
            if old_selearray is not None:
                sele_state[old_selearray] = selenr
            
        
    def PRIVATE_get_sele(self, sele, none_is_all=False):
        try:
            selenr = self.PRIVATE_sele.index(sele) + 1
        except ValueError:
            if sele == "sele":
                if none_is_all:
                    return slice(None, None)
                return None
            else:
                raise ValueError("Unknown selection '%s'" % sele)
        atomstate = self.atomstate.data
        selearray = (atomstate["sele"] == selenr)        
        return selearray

    def get_selection(self, sele="sele"):
        import numpy as np
        import pandas as pd 
        from collections import OrderedDict
        sele_array = self.PRIVATE_get_sele(sele)
        if sele_array is None:
            return
        atomstate = self.atomstate.data
        arr_atomstate = OrderedDict()
        for key in self.atomstate.dtype.fields.keys():
            if key in ("sele", "repr"):
                continue            
            v = atomstate[key]
            if key == "obj":
                names = [""] + self.PRIVATE_objects.data
                v = np.array(names)[v]
            elif v.dtype.kind == "S":
                v = v.astype("U")
            arr_atomstate[key] = v
        return pd.DataFrame(arr_atomstate).iloc[sele_array]

    def PRIVATE_change_repr(self, representation, sele, op):
        sele_array = self.PRIVATE_get_sele(sele, none_is_all=True)
        rep = self.atomstate.data["repr"]
        if representation == "all":
            assert op == "hide"
            rep[sele_array] = 0
            return
        reprnr = self.PRIVATE_repr.index(representation)
        repr_bit = 1 << reprnr
        if op == "show":
            rep[sele_array] |= repr_bit
        elif op == "hide":
            rep[sele_array] ^= repr_bit
        elif op == "show-as":
            rep[sele_array] = repr_bit

    def show(self, representation, sele="sele"):
        return self.PRIVATE_change_repr(representation, sele, "show")

    def hide(self, representation="all", sele="sele"):
        return self.PRIVATE_change_repr(representation, sele, "hide")

    def show_as(self, representation, sele="sele"):
        return self.PRIVATE_change_repr(representation, sele, "show-as")

    def ngl_representations(self):
        repr_mapping = {
            "cartoon": "cartoon", 
            "spacefill": "spacefill",
            "lines": "line", 
            "sticks": "licorice",
            "ball+sticks": "ball+stick",
            "dots": "point"  
        }
        import numpy as np
        result = []
        repr_state = self.atomstate.data["repr"]
        index_state = self.atomstate.data["index"]
        for n, rep in enumerate(self.PRIVATE_repr):
            reprbit = 1 << n
            indices = index_state[np.nonzero(repr_state & reprbit)]
            if not len(indices):
                continue
            result.append({
                "type": repr_mapping[rep],
                "params": {
                    "sele": "@" + ",".join([str(i) for i in indices])
                }
            })
        return result

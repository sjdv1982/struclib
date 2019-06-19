from seamless.silk import Silk
from seamless.silk.meta import meta, validator

class StructureState(metaclass=meta):
    
    version = "OK"

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
    # up to 255 objects   
    atomic_state_dtype =  \
        [("obj", 'uint8')] + \
        atomic_dtype + \
        [
            ("sele", 'uint64'),   # up to 64 selections
            ("repr", 'uint16'),    # up to 16 representations
            ("color", 'uint8')    # up to 256 colors
        ]

    # Names of colors and color schemes:
    PRIVATE_colors = [
        "atomindex",
        "bfactor",
        "chainid",
        "electrostatic",
        "element",
        "hydrophobicity",
        "modelindex",
        "occupancy",
        "residueindex",
        "resname",
        "sstruc"
    ] + \
    list({
        "aliceblue": "#f0f8ff",
        "antiquewhite": "#faebd7",
        "aqua": "#00ffff",
        "aquamarine": "#7fffd4",
        "azure": "#f0ffff",
        "beige": "#f5f5dc",
        "bisque": "#ffe4c4",
        "black": "#000000",
        "blanchedalmond": "#ffebcd",
        "blue": "#0000ff",
        "blueviolet": "#8a2be2",
        "brown": "#a52a2a",
        "burlywood": "#deb887",
        "cadetblue": "#5f9ea0",
        "chartreuse": "#7fff00",
        "chocolate": "#d2691e",
        "coral": "#ff7f50",
        "cornflowerblue": "#6495ed",
        "cornsilk": "#fff8dc",
        "crimson": "#dc143c",
        "cyan": "#00ffff",
        "darkblue": "#00008b",
        "darkcyan": "#008b8b",
        "darkgoldenrod": "#b8860b",
        "darkgray": "#a9a9a9",
        "darkgreen": "#006400",
        "darkgrey": "#a9a9a9",
        "darkkhaki": "#bdb76b",
        "darkmagenta": "#8b008b",
        "darkolivegreen": "#556b2f",
        "darkorange": "#ff8c00",
        "darkorchid": "#9932cc",
        "darkred": "#8b0000",
        "darksalmon": "#e9967a",
        "darkseagreen": "#8fbc8f",
        "darkslateblue": "#483d8b",
        "darkslategray": "#2f4f4f",
        "darkslategrey": "#2f4f4f",
        "darkturquoise": "#00ced1",
        "darkviolet": "#9400d3",
        "deeppink": "#ff1493",
        "deepskyblue": "#00bfff",
        "dimgray": "#696969",
        "dimgrey": "#696969",
        "dodgerblue": "#1e90ff",
        "firebrick": "#b22222",
        "floralwhite": "#fffaf0",
        "forestgreen": "#228b22",
        "fuchsia": "#ff00ff",
        "gainsboro": "#dcdcdc",
        "ghostwhite": "#f8f8ff",
        "gold": "#ffd700",
        "goldenrod": "#daa520",
        "gray": "#808080",
        "green": "#008000",
        "greenyellow": "#adff2f",
        "grey": "#808080",
        "honeydew": "#f0fff0",
        "hotpink": "#ff69b4",
        "indianred": "#cd5c5c",
        "indigo": "#4b0082",
        "ivory": "#fffff0",
        "khaki": "#f0e68c",
        "lavender": "#e6e6fa",
        "lavenderblush": "#fff0f5",
        "lawngreen": "#7cfc00",
        "lemonchiffon": "#fffacd",
        "lightblue": "#add8e6",
        "lightcoral": "#f08080",
        "lightcyan": "#e0ffff",
        "lightgoldenrodyellow": "#fafad2",
        "lightgray": "#d3d3d3",
        "lightgreen": "#90ee90",
        "lightgrey": "#d3d3d3",
        "lightpink": "#ffb6c1",
        "lightsalmon": "#ffa07a",
        "lightseagreen": "#20b2aa",
        "lightskyblue": "#87cefa",
        "lightslategray": "#778899",
        "lightslategrey": "#778899",
        "lightsteelblue": "#b0c4de",
        "lightyellow": "#ffffe0",
        "lime": "#00ff00",
        "limegreen": "#32cd32",
        "linen": "#faf0e6",
        "magenta": "#ff00ff",
        "maroon": "#800000",
        "mediumaquamarine": "#66cdaa",
        "mediumblue": "#0000cd",
        "mediumorchid": "#ba55d3",
        "mediumpurple": "#9370db",
        "mediumseagreen": "#3cb371",
        "mediumslateblue": "#7b68ee",
        "mediumspringgreen": "#00fa9a",
        "mediumturquoise": "#48d1cc",
        "mediumvioletred": "#c71585",
        "midnightblue": "#191970",
        "mintcream": "#f5fffa",
        "mistyrose": "#ffe4e1",
        "moccasin": "#ffe4b5",
        "navajowhite": "#ffdead",
        "navy": "#000080",
        "oldlace": "#fdf5e6",
        "olive": "#808000",
        "olivedrab": "#6b8e23",
        "orange": "#ffa500",
        "orangered": "#ff4500",
        "orchid": "#da70d6",
        "palegoldenrod": "#eee8aa",
        "palegreen": "#98fb98",
        "paleturquoise": "#afeeee",
        "palevioletred": "#db7093",
        "papayawhip": "#ffefd5",
        "peachpuff": "#ffdab9",
        "peru": "#cd853f",
        "pink": "#ffc0cb",
        "plum": "#dda0dd",
        "powderblue": "#b0e0e6",
        "purple": "#800080",
        "rebeccapurple": "#663399",
        "red": "#ff0000",
        "rosybrown": "#bc8f8f",
        "royalblue": "#4169e1",
        "saddlebrown": "#8b4513",
        "salmon": "#fa8072",
        "sandybrown": "#f4a460",
        "seagreen": "#2e8b57",
        "seashell": "#fff5ee",
        "sienna": "#a0522d",
        "silver": "#c0c0c0",
        "skyblue": "#87ceeb",
        "slateblue": "#6a5acd",
        "slategray": "#708090",
        "slategrey": "#708090",
        "snow": "#fffafa",
        "springgreen": "#00ff7f",
        "steelblue": "#4682b4",
        "tan": "#d2b48c",
        "teal": "#008080",
        "thistle": "#d8bfd8",
        "tomato": "#ff6347",
        "turquoise": "#40e0d0",
        "violet": "#ee82ee",
        "wheat": "#f5deb3",
        "white": "#ffffff",
        "whitesmoke": "#f5f5f5",
        "yellow": "#ffff00",
        "yellowgreen": "#9acd32"
    }.keys())

    def __init__(self):
        import numpy as np
        atomic_state_dtype = [tuple(x) for x in self.atomic_state_dtype]
        atomic_state_dtype = np.dtype(atomic_state_dtype, align=True)

        # Array of the atomic state
        self.atomstate = np.zeros(0, dtype=atomic_state_dtype)
        
        # Names of molecular objects
        self.PRIVATE_objects = []  
        
        # Names of representations
        self.PRIVATE_repr = [
            "cartoon", "ribbon", "rope",
            "base", "surface", "trace", "tube",
            "spacefill", "point", "line",
            "licorice", "ball+sticks", "hyperball",
            "label", "validation"
        ]  
        
        # Names of selections
        self.PRIVATE_sele = ["sele"]
        self.PRIVATE_active_selection = "sele"

    
    def parse_pdb(self, pdb, obj=None):
        import warnings
        import Bio.PDB
        from Bio.PDB.StructureBuilder import PDBConstructionWarning
        warnings.simplefilter('ignore', PDBConstructionWarning)
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
        a["color"] = self.PRIVATE_colors.index("chainid")
        a["obj"] = len(self.PRIVATE_objects)
        count = 0
        for modelnr, model in enumerate(struc.get_models()):
            atomlist = list(model.get_atoms())
            atomlist.sort(key=lambda atom: atom.serial_number)
            for atom in atomlist:
                residue = atom.get_parent()
                hetero, resid, icode = residue.get_id()
                segid = residue.segid
                resname = residue.resname
                chainid = residue.get_parent().id
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
        
    def select(self, query, sele=None):
        import seamless
        import numpy as np
        #Bizarre bug: the following code gives an error, but not if it is launched directly from the command line.
        # Must be in pandas/numexpr C code.
        # Hunches of possible causes: 1. something related to multiprocessing, 2. Something related to the bool datatype,
        #  3. Something related to globals() (local_dict solves the issue), 4. Something related to strides/memory alignment
        # For now, use local_dict and cast all bools to ints
        #
        #  selearray = seamless.pandeval.eval("sele", global_dict={"sele": np.array([0,1,1,0,0,1],bool)})
        if sele is None:
            sele = self.PRIVATE_active_selection
        self.PRIVATE_active_selection = sele                
        old_selearray = None
        sele_state = self.atomstate.data["sele"]
        try:
            selenr = self.PRIVATE_sele.index(sele) + 1                        
            selebit = np.uint64(1 << selenr)
            old_selearray = (sele_state & selebit)            
        except ValueError:
            self.PRIVATE_sele.append(sele)
            selenr = len(self.PRIVATE_sele)
            maxsel = 8 * sele_state.itemsize
            if selenr >= maxsel:
                raise ValueError("Maximum number of selections %d reached" % maxsel)
        selebit = np.uint64(1 << selenr)
        dic = {key: self.atomstate.data[key] for key in self.atomstate.dtype.fields.keys() \
               if key not in ("obj", "repr", "sele")}
        objects = np.array([""] + self.PRIVATE_objects.data).astype("S")
        dic["all"] = np.ones(len(self.atomstate))
        dic["obj"] = objects[self.atomstate.data["obj"]]
        dic["backbone"] = (self.atomstate.data["name"] == b"CA") | \
                (self.atomstate.data["name"] == b"C") | \
                (self.atomstate.data["name"] == b"O") | \
                (self.atomstate.data["name"] == b"N")
        for xselenr0, xsele in enumerate(self.PRIVATE_sele):
            xselenr = xselenr0 + 1
            if xsele in dic or (xselenr == selenr and old_selearray is None):
                continue
            xselebit = np.uint64(1 << xselenr)
            xselearray = (sele_state & xselebit).astype(bool).astype(np.uint8)
            dic[xsele] = xselearray
        try:
            selearray = seamless.pandeval.eval(query, local_dict=dic, align_result=False, str_as_bytes=True)
            selearray = np.ascontiguousarray(selearray).astype(bool).copy() # another Heisenbug; something is wrong with memory!
            print("%d atoms selected" % (selearray>0).sum())
            self.PRIVATE_unselect_noshift(sele)
            sele_state[selearray] |=  selebit
        except:
            import traceback; traceback.print_exc()
            if old_selearray is not None:
                sele_state[old_selearray] |= selebit
                    
    def PRIVATE_unselect_noshift(self, sele):
        try:
            selenr = self.PRIVATE_sele.index(sele) + 1
        except ValueError:
            raise ValueError("Unknown selection '%s'" % sele)
        atomstate_sele = self.atomstate.data["sele"]
        selebit = 1 << selenr
        dt = atomstate_sele.dtype.type
        cmpl_sele = ~dt(selebit)
        atomstate_sele &= cmpl_sele
        
    def unselect(self, sele=None):
        if sele is None:
            sele = self.PRIVATE_active_selection
        if sele == "sele":
            return self.PRIVATE_unselect_noshift("sele")
        try:
            selenr = self.PRIVATE_sele.index(sele) + 1
        except ValueError:
            raise ValueError("Unknown selection '%s'" % sele)
        atomstate_sele = self.atomstate.data["sele"]
        dt = atomstate_sele.dtype.type
        bits_to_keep = dt(2**selenr-1)
        bits_to_shift = ~bits_to_keep
        atomstate_sele[:] = ((atomstate_sele  >> 1) & bits_to_shift) | (atomstate_sele & bits_to_keep)
        self.PRIVATE_sele.remove(sele)
        if self.PRIVATE_active_selection == sele:
            self.PRIVATE_active_selection = "sele"

    def PRIVATE_get_sele(self, sele):
        try:
            selenr = self.PRIVATE_sele.index(sele) + 1
        except ValueError:
            raise ValueError("Unknown selection '%s'" % sele)
        atomstate = self.atomstate.data
        selebit = 1 << selenr
        selearray = (atomstate["sele"] & selebit).astype(bool)
        return selearray

    def get_selection(self, sele=None, format="mask"):
        assert format in ("mask", "pandas")
        import numpy as np
        import pandas as pd 
        from collections import OrderedDict
        if sele is None:
            sele = self.PRIVATE_active_selection
        mask = self.PRIVATE_get_sele(sele)
        if format == "mask":
            return mask
        atomstate = self.atomstate.data
        arr_atomstate = OrderedDict()
        for key in self.atomstate.dtype.fields.keys():
            if key in ("sele", "repr", "color"):
                continue                        
            v = atomstate[key]
            if key == "obj":
                names = [""] + self.PRIVATE_objects.data
                v = np.array(names)[v]
            elif v.dtype.kind == "S":
                v = v.astype("U")
            arr_atomstate[key] = v
        return pd.DataFrame(arr_atomstate).iloc[mask]

    def PRIVATE_change_repr(self, representation, sele, op):
        repr_mapping = {
            "lines": "line", 
            "sticks": "licorice",
            "ball+sticks": "ball+stick",
            "dots": "point"  
        }
        if representation in repr_mapping:
            representation = repr_mapping[representation]
        sele_array = self.PRIVATE_get_sele(sele)
        rep = self.atomstate.data["repr"]
        if representation == "all":
            assert op == "hide"
            rep[:] = 0
            return
        reprnr = self.PRIVATE_repr.index(representation)
        repr_bit = 1 << reprnr
        if op == "show":
            rep[sele_array] |= repr_bit
        elif op == "hide":
            rep[sele_array] ^= repr_bit
        elif op == "show-as":
            rep[sele_array] = repr_bit

    def show(self, representation, sele=None):
        if sele is None:
            sele = self.PRIVATE_active_selection
        return self.PRIVATE_change_repr(representation, sele, "show")

    def hide(self, representation="all", sele=None):
        if sele is None:
            sele = self.PRIVATE_active_selection
        return self.PRIVATE_change_repr(representation, sele, "hide")

    def show_as(self, representation, sele=None):
        if sele is None:
            sele = self.PRIVATE_active_selection
        return self.PRIVATE_change_repr(representation, sele, "show-as")

    def color(self, color, sele=None):
        if sele is None:
            sele = self.PRIVATE_active_selection
        if color not in self.PRIVATE_colors:
            raise ValueError("Unknown color '%s'" % color)
        color_index = self.PRIVATE_colors.index(color)
        sele_array = self.PRIVATE_get_sele(sele)
        self.atomstate.data["color"][sele_array] = color_index

    def ngl_representations(self):
        import numpy as np
        results = {}
        for objnr, obj in enumerate(self.PRIVATE_objects):
            result = []            
            repr_state = self.atomstate.data["repr"]
            color_state = self.atomstate.data["color"]
            obj_mask = (self.atomstate.data["obj"] == objnr + 1)
            repr_state = repr_state[obj_mask]
            color_state = color_state[obj_mask]
            colors = np.unique(color_state)
            for n, rep in enumerate(self.PRIVATE_repr):
                reprbit = 1 << n
                repr_mask = (repr_state & reprbit).astype(bool)
                if not repr_mask.sum():
                    continue
                for color in colors:
                    color_mask = (color_state == color)
                    indices = np.nonzero(color_mask & repr_mask)[0]
                    if not len(indices):
                        continue
                    color_name = self.PRIVATE_colors[color]
                    result.append({
                        "type": rep,
                        "params": {
                            "sele": "@" + ",".join([str(i) for i in indices]),
                            "color": color_name
                        }
                    })
            if len(result):
                results[obj] = result
        return results

test = StructureState() #for type inference

result = StructureState.schema.dict # we have to put this file in a Transformer, until Seamless will support high-level Python modules

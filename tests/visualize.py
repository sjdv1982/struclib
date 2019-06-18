from seamless.silk import Silk
import requests
StructureState = Silk(schema=structurestate_schema)
struc = StructureState()
struc.data = struc_data

struc.hide()
struc.select('obj == "1ACB" and chain == "E"')
struc.show_as("hyperball")
struc.color("orange")
struc.select('obj == "1ACB" and chain == "I"')
struc.show_as("cartoon")
struc.color("green")
struc.select('obj == "1ACB" and chain == "I" and resid > 40 and resid < 56', 'activesite')
struc.color("blue")
struc.show("licorice")
struc.show("label")
result = struc.ngl_representations()
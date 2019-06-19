from seamless.silk import Silk
StructureState = Silk(schema=structurestate_schema)
struc = StructureState()
struc.data = struc_data

for lineno0, line in enumerate(visualization.splitlines()):
    lineno = lineno0 + 1
    pound = line.find("#")
    if pound > -1:
        line = line[:pound]
    if not len(line.strip()):
        continue
    terms = [t.strip() for t in line.split(",")]
    try:
        firstargs = terms[0].split()
        command = firstargs[0]
        args = []
        if len(firstargs) > 1:            
            args.append(" ".join(firstargs[1:]))
        args += terms[1:]
    except:
        raise SyntaxError("Line %d: %s" % (lineno, line))
    try:
        func = getattr(struc, command)
    except AttributeError:
        raise SyntaxError("Line %d, unknown command '%s'" % (lineno, command))
    func(*args)
    
result = {
    "mask": struc.get_selection(format="mask"),
    "table": struc.get_selection(format="pandas").to_html(),
    "ngl_representations": struc.ngl_representations(),
}
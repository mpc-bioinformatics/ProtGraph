def get_content(text, beginning="(", delimiter="->"):
    # A Sequence (or AA) is added
    # Get   X -> Y   Information
    # TODO duplicated in VAR_SEQ?
    idx = text.find(beginning)
    if idx != -1:
        text = text[:idx]
    xy = text.split(delimiter, 2)
    y_s = xy[1].strip().replace(" ", "")
    return y_s

def get_isoforms(text):
    isoforms = text[
        text.index("(") + 1 + 3: text.rfind(")")  # +3 to remove "in "
    ]
    isoforms = isoforms.replace("isoform ", "").replace(" and ", ", ")
    output_isoforms = []
    for isoform in isoforms.split(", "):
        isoform = isoform.strip()
        output_isoforms.append(", isoform " + isoform)
    return output_isoforms

def _get_qualifiers(edge):
    """ A simple method to retrieve qualifiers. It always returns a list """
    qualifiers = edge["qualifiers"] if "qualifiers" in edge.attributes() else []
    qualifiers = [] if qualifiers is None else qualifiers
    return qualifiers

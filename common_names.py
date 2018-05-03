from periodictable import formula
# lisst of common names with chemical formula

table = [("acac", "C5H7O2"),
         ("Et", "C2H5"),
         ("iPr", "C3H7"),
         ("TBAB", "C16H36BrN"),
         ('CTAB', 'C19H42BrN'),
         ("hz", "N2H4"),
         ("pKPt", "K2PtCl6"),
         ("pHPt", "H2PtCl6+5H2O")]

table = dict(table)


def parse_common_name(s):
    res = s
    for k, v in table.items():
        if k in s:
            res = s.replace(k, v)
    return res


def parse_formula(name):
    return formula(parse_common_name(name))

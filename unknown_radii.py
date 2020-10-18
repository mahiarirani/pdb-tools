def get_data(atom):
    if atom.element == "H":
        atom.radius = 1.0
        atom.polar = False

    elif atom.element == "C":
        atom.radius = 1.8
        atom.polar = False

    elif atom.element == "N":
        atom.radius = 1.6
        atom.polar = True

    elif atom.element == "O":
        atom.radius = 1.4
        atom.polar = True

    elif atom.element == "P":
        atom.radius = 1.9
        atom.polar = True

    elif atom.element == "S":
        atom.radius = 1.85
        atom.polar = False

    elif atom.element == "SE":
        atom.radius = 1.9
        atom.polar = False

    elif atom.element == "I":
        atom.radius = 2.094
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "F":
        atom.radius = 1.560
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Br" or atom.element == "BR":
        atom.radius = 1.978
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Cl" or atom.element == "CL":
        atom.radius = 1.74
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Al" or atom.element == "AL":
        atom.radius = 0.54
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "As" or atom.element == "AS":
        atom.radius = 0.58
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Au" or atom.element == "AU":
        atom.radius = 1.37
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Ba" or atom.element == "BA":
        atom.radius = 1.35
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Be" or atom.element == "BE":
        atom.radius = 0.45
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Bi" or atom.element == "BI":
        atom.radius = 1.03
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Ca" or atom.element == "CA":
        atom.radius = 1.00
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Cd" or atom.element == "CD":
        atom.radius = 0.95
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Cr" or atom.element == "CR":
        atom.radius = 0.73
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Cs" or atom.element == "CS":
        atom.radius = 1.67
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Cu" or atom.element == "CU":
        atom.radius = 0.73
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Fe" or atom.element == "FE":
        atom.radius = 0.65
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Ga" or atom.element == "GA":
        atom.radius = 0.62
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Ge" or atom.element == "GE":
        atom.radius = 0.73
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Hg" or atom.element == "HG":
        atom.radius = 1.02
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "K" or atom.element == "K":
        atom.radius = 1.38
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Li" or atom.element == "LI":
        atom.radius = 0.76
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Mg" or atom.element == "MG":
        atom.radius = 0.72
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Mn" or atom.element == "MN":
        atom.radius = 0.83
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Mo" or atom.element == "MO":
        atom.radius = 0.69
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Na" or atom.element == "NA":
        atom.radius = 1.02
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Ni" or atom.element == "NI":
        atom.radius = 0.69
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Pb" or atom.element == "PB":
        atom.radius = 1.19
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Pd" or atom.element == "PD":
        atom.radius = 0.86
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Rb" or atom.element == "RB":
        atom.radius = 1.52
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Sb" or atom.element == "SB":
        atom.radius = 0.76
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Sc" or atom.element == "SC":
        atom.radius = 0.75
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Sn" or atom.element == "SN":
        atom.radius = 0.69
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Sr" or atom.element == "SR":
        atom.radius = 1.18
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Tc" or atom.element == "TC":
        atom.radius = 0.65
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Ti" or atom.element == "TI":
        atom.radius = 0.86
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "V":
        atom.radius = 0.79
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Zn" or atom.element == "ZN":
        atom.radius = 0.74
        atom.polar = True if atom.name == atom.element else False

    elif atom.element == "Zr" or atom.element == "ZR":
        atom.radius = 0.72
        atom.polar = True if atom.name == atom.element else False

    else:
        atom.radius = 0
        atom.polar = False

    return atom

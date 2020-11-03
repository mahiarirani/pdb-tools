def get_data(element, name):
    if element == "H":
        radius = 1.0
        polar = False

    elif element == "C":
        radius = 1.8
        polar = False

    elif element == "N":
        radius = 1.6
        polar = True

    elif element == "O":
        radius = 1.4
        polar = True

    elif element == "P":
        radius = 1.9
        polar = True

    elif element == "S":
        radius = 1.85
        polar = False

    elif element == "SE":
        radius = 1.9
        polar = False

    elif element == "I":
        radius = 2.094
        polar = True if name == element else False

    elif element == "F":
        radius = 1.560
        polar = True if name == element else False

    elif element == "Br" or element == "BR":
        radius = 1.978
        polar = True if name == element else False

    elif element == "Cl" or element == "CL":
        radius = 1.74
        polar = True if name == element else False

    elif element == "Al" or element == "AL":
        radius = 0.54
        polar = True if name == element else False

    elif element == "As" or element == "AS":
        radius = 0.58
        polar = True if name == element else False

    elif element == "Au" or element == "AU":
        radius = 1.37
        polar = True if name == element else False

    elif element == "Ba" or element == "BA":
        radius = 1.35
        polar = True if name == element else False

    elif element == "Be" or element == "BE":
        radius = 0.45
        polar = True if name == element else False

    elif element == "Bi" or element == "BI":
        radius = 1.03
        polar = True if name == element else False

    elif element == "Ca" or element == "CA":
        radius = 1.00
        polar = True if name == element else False

    elif element == "Cd" or element == "CD":
        radius = 0.95
        polar = True if name == element else False

    elif element == "Cr" or element == "CR":
        radius = 0.73
        polar = True if name == element else False

    elif element == "Cs" or element == "CS":
        radius = 1.67
        polar = True if name == element else False

    elif element == "Cu" or element == "CU":
        radius = 0.73
        polar = True if name == element else False

    elif element == "Fe" or element == "FE":
        radius = 0.65
        polar = True if name == element else False

    elif element == "Ga" or element == "GA":
        radius = 0.62
        polar = True if name == element else False

    elif element == "Ge" or element == "GE":
        radius = 0.73
        polar = True if name == element else False

    elif element == "Hg" or element == "HG":
        radius = 1.02
        polar = True if name == element else False

    elif element == "K" or element == "K":
        radius = 1.38
        polar = True if name == element else False

    elif element == "Li" or element == "LI":
        radius = 0.76
        polar = True if name == element else False

    elif element == "Mg" or element == "MG":
        radius = 0.72
        polar = True if name == element else False

    elif element == "Mn" or element == "MN":
        radius = 0.83
        polar = True if name == element else False

    elif element == "Mo" or element == "MO":
        radius = 0.69
        polar = True if name == element else False

    elif element == "Na" or element == "NA":
        radius = 1.02
        polar = True if name == element else False

    elif element == "Ni" or element == "NI":
        radius = 0.69
        polar = True if name == element else False

    elif element == "Pb" or element == "PB":
        radius = 1.19
        polar = True if name == element else False

    elif element == "Pd" or element == "PD":
        radius = 0.86
        polar = True if name == element else False

    elif element == "Rb" or element == "RB":
        radius = 1.52
        polar = True if name == element else False

    elif element == "Sb" or element == "SB":
        radius = 0.76
        polar = True if name == element else False

    elif element == "Sc" or element == "SC":
        radius = 0.75
        polar = True if name == element else False

    elif element == "Sn" or element == "SN":
        radius = 0.69
        polar = True if name == element else False

    elif element == "Sr" or element == "SR":
        radius = 1.18
        polar = True if name == element else False

    elif element == "Tc" or element == "TC":
        radius = 0.65
        polar = True if name == element else False

    elif element == "Ti" or element == "TI":
        radius = 0.86
        polar = True if name == element else False

    elif element == "V":
        radius = 0.79
        polar = True if name == element else False

    elif element == "Zn" or element == "ZN":
        radius = 0.74
        polar = True if name == element else False

    elif element == "Zr" or element == "ZR":
        radius = 0.72
        polar = True if name == element else False

    else:
        radius = 0
        polar = False

    return radius, polar

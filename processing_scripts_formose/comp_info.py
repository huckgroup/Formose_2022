from pathlib import Path
import numpy as np
import matplotlib as mpl


class information:
    def __init__(self, comp_info):
        def colorFader(
            c1, c2, mix=0
        ):  # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
            c1 = np.array(mpl.colors.to_rgb(c1))
            c2 = np.array(mpl.colors.to_rgb(c2))
            return mpl.colors.to_hex((1 - mix) * c1 + mix * c2)

        path = Path(comp_info)
        path_file = path / "compound_information.csv"
        with open(path_file, "r") as f:
            f_n = []
            for line in f:
                f_n.append(line.strip("\n").split(","))

            ind = []
            working_name = []
            names = []
            SMILES = []
            colour = []
            for a, x in enumerate(f_n):
                if a > 0:
                    ind.append(x[0])
                    working_name.append(x[1])
                    name = x[2]
                    name = name.replace("alpha-beta", "\u03B1,\u03B2")
                    name = name.replace("alpha", "\u03B1")
                    name = name.replace("beta", "\u03B2")
                    name = name.replace("RRR", "R,R,R")
                    name = name.replace("RRS", "R,R,S")
                    name = name.replace("RSR", "R,S,R")
                    name = name.replace("RSS", "R,S,S")
                    name = name.replace("RR", "R,R")
                    name = name.replace("RS", "R,S")
                    names.append(name)
                    SMILES.append(x[3])
                    c1 = x[4]
                    c2 = "#c5c9c7"
                    c = []
                    n = 3
                    for x in range(n + 1):
                        c.append(colorFader(c1, c2, x / n))
                    colour.append(c)

        self.ind = ind
        self.working_name = working_name
        self.name = names
        self.SMILES = SMILES
        self.colour = colour


def load_colours_dict(filename):
    """
    Load a dictionary of compound colours from a .csv file.

    Parameters
    ----------
    filename: str or pathlib.Path
        Path to file.
    index: int
        Column containing the colour.

    Returns
    -------
    colours: dict()
    """

    with open(filename, "r") as file:
        text = file.read()

    lines = [x for x in text.split("\n") if x != ""]

    colours = dict()
    for line in lines[1:]:
        entries = line.split(",")
        colours[entries[3]] = entries[4]

    return colours

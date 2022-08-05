"""
directory to modules and data folder
"""


def load_config(filename):
    """
    Load a config file.

    Parameters
    ----------
    filename: st
        Path to the config file.

    Returns
    -------
    config: dict
    """

    with open(filename, "r") as file:
        text = file.read()

    lines = [x for x in text.split("\n") if x != ""]

    config = dict()
    for line in lines:
        tokens = line.split(",")
        config[tokens[0]] = tokens[1]

    return config


def load_exp_info(paths_folder, experiments):
    with open(paths_folder, "r") as f:
        exp = []
        for line in f:
            line = list(line.strip("\n").split(","))
            if line[0] in experiments:
                exp.append(line)
        list_exp = {}
        condition = []
        for x in exp:
            name = x[4]
            name = name.replace("stdev", "\u03C3")
            list_exp[x[0]] = []
            list_exp[x[0]] += [x[1], x[2], x[3], name, x[5]]
            condition.append(x[0])
        return list_exp, condition

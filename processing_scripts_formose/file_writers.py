def write_average_stdev_csv(averages, st_devs, comp_ind, filename=""):
    """
    Write averages and standard deviation dicts to a .csv file.

    Parameters
    ----------
    averages: dict()
    st_devs: dict()
    filename: str or pathlib.Path

    Returns
    -------
    None
    """

    data_string = "index,compound,average/ M,standard deviation/ M\n"

    for compound in averages:
        compound_token = compound.split("/")[0]
        ind = str(comp_ind[compound_token])
        data_string += f"{ind},{compound_token},{averages[compound]},"
        data_string += f"{st_devs[compound]}\n"

    with open(filename, "w") as file:
        file.write(data_string)
    print("Results written to output file: ", f"{filename}")


def write_corr_csv(corr, t_interval, ind, filename=""):

    data_string = "compound_ind,"
    for n in ind:
        data_string += f"{n},"
    data_string += f"\n"

    for a, x in enumerate(corr):
        data_string += f"{t_interval[a]}_s,"
        for y in x:
            data_string += f"{y[0]},"
        data_string += f"\n"

    for a, x in enumerate(corr):
        data_string += f"{t_interval[a]}_s_hex_col,"
        for y in x:
            data_string += f"{y[1]},"
        data_string += f"\n"

    with open(filename, "w") as file:
        file.write(data_string)
    print("Results written to output file: ", f"{filename}")


def write_rel_diff_csv(dic_rel_diff, exp, filename=""):
    """
    Write relative difference in concentration dicts to a .csv file.

    Parameters
    ----------
    dic_rel_diff: dict()
    filename: str or pathlib.Path

    Returns
    -------
    None
    """
    data_string = "compound,"

    for x in exp:
        data_string += f"{x},"
    data_string += f"\n"

    for compound in dic_rel_diff:
        data_string += f"{compound},"
        for x in dic_rel_diff[compound]:
            data_string += f"{x},"
        data_string += f"\n"

    with open(filename, "w") as file:
        file.write(data_string)
    print("Results written to output file: ", f"{filename}")

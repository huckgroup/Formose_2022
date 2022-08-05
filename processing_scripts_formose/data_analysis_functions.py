import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.colors as mcolors


def data_averages(data_report):
    """
    Calculate the averages of the compound traces in a data report object.

    Parameters
    ----------
    data_report: data_report.data_report

    Returns
    -------
    averages: dict()
    """

    averages = dict()

    for compound in data_report.data:
        averages[compound] = np.average(data_report.data[compound])

    return averages


def data_standard_deviations(data_report):
    """
    Calculate the standard deviations of the compound traces in a data report object.
    Parameters
    ----------
    data_report: data_report.data_report

    Returns
    -------
    st_devs: dict()
    """

    st_devs = dict()

    for compound in data_report.data:
        st_devs[compound] = np.std(data_report.data[compound], ddof=1)

    return st_devs


def data_p_values(data_report_1, data_report_2, list_comp):
    """
    Calculate the p-values for the compound traces in a data report.

    Parameters
    ----------
    data_report_1: data_report.data_report
    data_report_2: data_report.data_report

    Returns
    -------
    ttest: dict()
    """
    p_values = dict()

    steady_state_comp = [*data_report_1.data]
    perturbed_state_comp = [*data_report_2.data]

    steady_state = data_report_1.data
    perturbed_state = data_report_2.data

    dist_1 = []
    dist_2 = []

    # if a compound is not present in a data report, its series is defined to
    # be a sequence of zeros
    for _, compound in list_comp:
        token = "no_comp"
        # find steady state concentration for compound, zero values are removed from list
        for x in steady_state_comp:
            if compound in x:
                token = x
        if token in steady_state_comp:
            dist_1 = steady_state[token]
        else:
            dist_1 = [0]

        # find perturbed state concentration for compound, zero values are removed from list
        token = "no_comp"
        for x in perturbed_state_comp:
            if compound in x:
                token = x
        if token in perturbed_state_comp:
            dist_2 = perturbed_state[token]
        else:
            dist_2 = [0]

        result = stats.ttest_ind(dist_1, dist_2, nan_policy="omit")

        p_values[token] = result.pvalue
    return p_values


def differential_means(data, t_interval, sample_time, l):
    """
    Find the variation in the signal on timescale t_interval.

    Parameters
    ----------
    data: data_report.data_report
    t_interval: 1D list
    sample_time: 1D list
    l = 1D list

    Returns
    -------
    differentials: list of list of lists with respectively time intervals, compound differential and values in the differential
    """
    compounds = [*data]

    # def differential_means(val,t_interval,sample_time):
    differentials = []
    for x in t_interval:
        comb_lag = []
        interval = x / sample_time
        start = interval - 1

        token = ""
        for _, y in l:
            condition = False
            for key in compounds:
                if y in key:
                    token = key
                    condition = True
            if condition == True:
                d = data[token]
                lag = []
                for a, z in enumerate(d):
                    if a > start:
                        sample_range = d[int(a - interval) : int(a)]
                        lag.append(np.mean(sample_range))

                d_lag = []
                for a, z in enumerate(lag):
                    if a > start:
                        d_lag.append(z - lag[int(a - interval)])

                m_lag = np.mean(d_lag)

                zero_a = []
                for x in d_lag:
                    zero_a.append(x - m_lag)

                norm_flow = zero_a / np.std(zero_a)

                comb_lag.append(norm_flow)
        differentials.append(comb_lag)
    return differentials


def correlation(val, flow, interval):
    """
    Find the Pearson correlation between the differential of the Ca(OH)2 input flow and compound output
    at different time intervals.

    Parameters
    ----------
    val: list of list of lists with respectively time intervals, compound differential and values in the differential
    flow: list of list of lists with respectively time intervals, input flow differential and values in the differential
    interval: list with time intervals

    Returns
    -------
    correlation: list of list of list with respectively time interval and correlation + hex color per compound to the input
    """

    ### calculate Pearson correlation at different time intervals betweeen input flow and output compounds ###
    all_corr = []
    ### iterates through differential at time scales of time interval ###
    for a, x in enumerate(flow):
        corr = []
        ### iterates through compounds for respective time scale of flow ###
        for y in val[a]:
            ### Pearson correlation ###
            c = stats.pearsonr(y, x[0])
            corr.append(c[0])
        all_corr.append(corr)

    ### calculate largest positive/negative correlation value in data for color scale ###
    corr_range = []
    flat = np.array(all_corr).flatten()
    if np.amax(flat) > abs(np.amin(flat)):
        value = np.amax(flat)
        corr_range = [-value, value]
    elif abs(np.amin(flat)) > np.amax(flat):
        value = abs(np.amin(flat))
        corr_range = [-value, value]

    ### color map to express calculated correlation value ###
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "",
        ["midnightblue", "#0000D6", "lightskyblue", "white", "pink", "red", "maroon"],
        N=360,
    )

    ### multiplication factor for to  calculate color for each correlation value from 'cmap' ###
    mult = 360 / (abs(corr_range[0]) + corr_range[1])

    ### calculate colours for each correlation value and put each color in list with corresponding correlation value ###
    correlation = []
    for a, x in enumerate(interval):
        ### list of correlations for each time interval ###
        l_col = []
        for y in all_corr[a]:
            ### calculate color value from correlation value ###
            conv_val = 360 - (y + corr_range[1]) * mult
            ### create list with correlation value and corresponding hex color ###
            l_col.append([y, mcolors.rgb2hex(cmap(360 - int(conv_val)))])
        ### list that contains lists for different time scales ###
        correlation.append(l_col)

    return correlation


def difference_average(data_report_1, data_report_2, list_comp, dic_rel_diff):
    """
    Calculate relative difference between perturbed state versus the steady state compound average in a data report.

    Parameters
    ----------
    data_report_1: data_report.data_report
    data_report_2: data_report.data_report
    list_comp: list
    dic_rel_diff: dict of lists

    Returns
    -------
    relative difference per compound: dict()
    """
    steady_state_comp = [*data_report_1.data]
    perturbed_state_comp = [*data_report_2.data]

    steady_state = data_report_1.data
    perturbed_state = data_report_2.data

    mean_1 = []
    mean_2 = []

    # if a compound is not present in a data report, its series is defined to
    # be a sequence of zeros
    for ind, compound in list_comp:
        token = "no_comp"
        # find steady state concentration for compound, zero values are removed from list
        for x in steady_state_comp:
            if compound in x:
                token = x
        if token in steady_state_comp:
            mean_1 = steady_state[token]
            mean_1 = [i for i in mean_1 if i != 0]
        else:
            mean_1 = [0]

        # find perturbed state concentration for compound, zero values are removed from list
        token = "no_comp"
        for x in perturbed_state_comp:
            if compound in x:
                token = x
        if token in perturbed_state_comp:
            mean_2 = perturbed_state[token]
            mean_2 = [i for i in mean_2 if i != 0]
        else:
            mean_2 = np.nan

        # calculate relative difference from with steady state 'mean_1' and perturbed state 'mean_2'
        result = (np.mean(mean_2) - np.mean(mean_1)) / np.mean(mean_1)
        if np.isnan(result) == True:
            result = 0

        dic_rel_diff[ind].append(result)

    return dic_rel_diff


def normalized_difference(dic_rel_diff, list_comp):
    """
    Normalizes the relative difference from steady steate for each of the observed compounds over EXP001 - EXP012

    Parameters
    ----------
    dic_rel_diff: data_report.data_report
    data_report_2: data_report.data_report

    Returns
    -------
    relative difference per compound: dict()
    """
    # calculate normalized difference for each compound
    for ind, _ in list_comp:
        c_list = dic_rel_diff[ind]
        num_list = []
        # exclude inf for finding largest difference from steady state
        for i in c_list:
            if np.isinf(i) == False:
                num_list.append(i)

        # calculate minumum and maximum value
        l_max = np.amax(num_list)
        l_min = np.amin(num_list)

        # set l_max to largest difference in list
        if abs(l_min) > l_max:
            l_max = abs(l_min)

        # normalize list for compound with l_max
        norm = []
        for i in c_list:
            if np.isinf(i) == True:
                norm.append(1)
            else:
                norm.append(i / l_max)
        # for compound replace relative difference list with normalized list in dictionary
        dic_rel_diff[ind] = norm
    return dic_rel_diff

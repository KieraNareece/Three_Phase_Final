# ======================================================================================================================
# Script Name: Function.py
# File Name: Three_Phase_Final
# Date: 06-09-2022
# Creator: K.N. Lambrecht
# Purpose: This script is used for storage of useful parameters
# ======================================================================================================================

# Import libraries to use in the following code

import types
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.optimize import least_squares
from matplotlib import pyplot as plt
import random


def vector_to_namespace(vector, var_list):
    # Vector:   A vector of parameters
    # Var_list: A list of parameter names

    # This function is used to convert the vectors used in the regression to a namespace
    # This store the vector into a holding dictionary which is then converted to a namespace
    # The dictionary elements are defined by the list of variables

    dictionary = {var_list[i]: vector[i] for i in range(len(var_list))}
    return types.SimpleNamespace(**dictionary)


def namespace_to_vector(namesp, var_list):
    # namesp:   A namespace of parameters
    # var_list: A list of parameter names

    # This function is used to convert a namespace into a vector for use in the regression
    # The namespace is first converted into a dictionary and then a vector is populated for this

    dictionary = namesp.__dict__
    vector = np.zeros(len(var_list))
    for i in range(len(var_list)):
        vector[i] = dictionary[var_list[i]]
    return vector


def rename_dataframe(df, index):
    # df:       A namespace of data frames
    # index:    A list containing data frame names
    # As the files store the experimental data with column names that have white space
    # these need to be renamed to be used in the code.

    dictionary = df.__dict__

    for i in range(0, len(index)):
        a = list(dictionary[index[i]].columns.values)
        rename_dic = {a[j]: (a[j] + ' ' + index[i]).replace(' ', '_') for j in range(len(a))}
        dictionary[index[i]].rename(columns=rename_dic, inplace=True)

    ns = types.SimpleNamespace(**dictionary)
    return ns


def experimental_data(ns, index):
    # ns:      A namespace of dataframes
    # index:   A list of data frame names
    # The files contain each repeat separately
    # For the purposes of the regression, the repeats need to be averaged which this function performs
    # In addition to this, the files containing the lees data refer to the composition of the must after lees contact,
    # and this means that the relevant days must be subtracted from the must values and multiplied by the volume to
    # obtain the total mass of a specific compound moving into the lees phase

    hold = rename_dataframe(ns, index)
    dictionary = hold.__dict__
    new_dictionary = {}
    values = {}

    for i in range(0, len(index)):
        a = list(dictionary[index[i]].columns.values)
        a.remove('Hour_' + index[i])
        values[index[i]] = a

        array = dictionary[index[i]].to_records(index=False)
        s = int(len(a))
        exp_df = {}

        for j in range(0, s):
            exp_df[a[j]] = array[a[j]]

        new_ns = types.SimpleNamespace(**exp_df)
        new_dictionary[index[i]] = new_ns

    data = types.SimpleNamespace(**new_dictionary)

    return data, values


def standard_deviations(df_namespace, phase_name, p):
    # df_namespace: A namespace of sorted dataframes
    # phase_name:   A list of dataframe names
    # p:            A namespace of parameters
    # The function uses lists contained in p, to call certain columns from the data frame and calculate the mean and
    # standard deviation of the experimental data

    df = df_namespace.__dict__
    phase = phase_name

    stdev = {}
    mean = {}

    for i in range(0, len(p.scenario)):

        if p.scenario[i] == 'Cultivar':
            for j in range(0, len(phase)):
                if phase[j] == 'Cap':
                    for k in range(0, len(p.cultivars)):
                        for l in range(0, len(p.cap_components)):
                            view_data = (df['Cap'].filter(regex=p.cultivars[k]).filter(regex=p.cap_components[l],
                                                                                       axis=1)).transpose()
                            view_data_std = (view_data.std()).to_numpy()
                            view_data_mean = (view_data.mean()).to_numpy()
                            stdev[(phase[j] + '_' + p.cultivars[k] + '_' + p.cap_components[l])] = view_data_std
                            mean[(phase[j] + '_' + p.cultivars[k] + '_' + p.cap_components[l])] = view_data_mean

                if phase[j] == 'Must':
                    for k in range(0, len(p.cultivars)):
                        for l in range(0, len(p.must_components)):
                            view_data = (df['Must'].filter(regex=p.cultivars[k]).filter(regex=p.must_components[l],
                                                                                        axis=1)).transpose()
                            view_data_std = (view_data.std()).to_numpy()
                            view_data_mean = (view_data.mean()).to_numpy()
                            stdev[(phase[j] + '_' + p.cultivars[k] + '_' + p.must_components[l])] = view_data_std
                            mean[(phase[j] + '_' + p.cultivars[k] + '_' + p.must_components[l])] = view_data_mean

                if phase[j] == 'Lees':
                    for k in range(0, len(p.cultivars)):
                        for l in range(0, len(p.lees_components)):
                            view_data = (df['Lees'].filter(regex=p.cultivars[k]).filter(regex=p.lees_components[l],
                                                                                        axis=1)).transpose()
                            view_data_std = (view_data.std()).to_numpy()
                            view_data_mean = (view_data.mean()).to_numpy()
                            stdev[(phase[j] + '_' + p.cultivars[k] + '_' + p.lees_components[l])] = view_data_std
                            mean[(phase[j] + '_' + p.cultivars[k] + '_' + p.lees_components[l])] = view_data_mean

        if p.scenario[i] == 'Temperature':
            for j in range(0, len(phase)):
                if phase[j] == 'Cap':
                    for k in range(0, len(p.temperature)):
                        for l in range(0, len(p.cap_components)):
                            view_data = (df['Cap'].filter(regex=p.temperature[k]).filter(regex=p.cap_components[l],
                                                                                         axis=1)).transpose()
                            view_data_std = (view_data.std()).to_numpy()
                            view_data_mean = (view_data.mean()).to_numpy()
                            stdev[(phase[j] + '_' + p.temperature[k] + '_' + p.cap_components[l])] = view_data_std
                            mean[(phase[j] + '_' + p.temperature[k] + '_' + p.cap_components[l])] = view_data_mean

                if phase[j] == 'Must':
                    for k in range(0, len(p.temperature)):
                        for l in range(0, len(p.must_components)):
                            view_data = (df['Must'].filter(regex=p.temperature[k]).filter(regex=p.must_components[l],
                                                                                          axis=1)).transpose()
                            view_data_std = (view_data.std()).to_numpy()
                            view_data_mean = (view_data.mean()).to_numpy()
                            stdev[(phase[j] + '_' + p.temperature[k] + '_' + p.must_components[l])] = view_data_std
                            mean[(phase[j] + '_' + p.temperature[k] + '_' + p.must_components[l])] = view_data_mean

                if phase[j] == 'Lees':
                    for k in range(0, len(p.temperature)):
                        for l in range(0, len(p.lees_components)):
                            view_data = (df['Lees'].filter(regex=p.temperature[k]).filter(regex=p.lees_components[l],
                                                                                          axis=1)).transpose()
                            view_data_std = (view_data.std()).to_numpy()
                            view_data_mean = (view_data.mean()).to_numpy()
                            stdev[(phase[j] + '_' + p.temperature[k] + '_' + p.lees_components[l])] = view_data_std
                            mean[(phase[j] + '_' + p.temperature[k] + '_' + p.lees_components[l])] = view_data_mean

        if p.scenario[i] == 'Specific':
            for j in range(0, len(phase)):
                if phase[j] == 'Cap':
                    for k in range(0, len(p.cultivars)):
                        for l in range(0, len(p.temperature)):
                            for m in range(0, len(p.cap_components)):
                                view_data = (df['Cap'].filter(regex=p.cultivars[k]).filter(
                                    regex=p.temperature[l]).filter(regex=p.cap_components[m], axis=1)).transpose()
                                view_data_std = (view_data.std()).to_numpy()
                                view_data_mean = (view_data.mean()).to_numpy()
                                stdev[(phase[j] + '_' + p.cultivars[k] + '_' + p.temperature[l] + '_' +
                                       p.cap_components[m])] = view_data_std
                                mean[(phase[j] + '_' + p.cultivars[k] + '_' + p.temperature[l] + '_' +
                                      p.cap_components[m])] = view_data_mean

                if phase[j] == 'Must':
                    for k in range(0, len(p.cultivars)):
                        for l in range(0, len(p.temperature)):
                            for m in range(0, len(p.must_components)):
                                view_data = (df['Must'].filter(regex=p.cultivars[k]).filter(
                                    regex=p.temperature[l]).filter(regex=p.must_components[m], axis=1)).transpose()
                                view_data_std = (view_data.std()).to_numpy()
                                view_data_mean = (view_data.mean()).to_numpy()
                                stdev[(phase[j] + '_' + p.cultivars[k] + '_' + p.temperature[l] + '_' +
                                       p.must_components[m])] = view_data_std
                                mean[(phase[j] + '_' + p.cultivars[k] + '_' + p.temperature[l] + '_' +
                                      p.must_components[m])] = view_data_mean

                if phase[j] == 'Lees':
                    for k in range(0, len(p.cultivars)):
                        for l in range(0, len(p.temperature)):
                            for m in range(0, len(p.lees_components)):
                                view_data = (df['Lees'].filter(regex=p.cultivars[k]).filter(
                                    regex=p.temperature[l]).filter(regex=p.lees_components[m], axis=1)).transpose()
                                view_data_std = (view_data.std()).to_numpy()
                                view_data_mean = (view_data.mean()).to_numpy()
                                stdev[(phase[j] + '_' + p.cultivars[k] + '_' + p.temperature[l] + '_' +
                                       p.lees_components[m])] = view_data_std
                                mean[(phase[j] + '_' + p.cultivars[k] + '_' + p.temperature[l] + '_' +
                                      p.lees_components[m])] = view_data_mean

    stdev_hold = types.SimpleNamespace(**stdev)
    mean_hold = types.SimpleNamespace(**mean)

    return stdev_hold, mean_hold


def int_var(t, x, p):
    # t - Time span of fermentation
    # x - A name space of state variables used in conjunction with the stoichiometric matrix
    # p - A namespace containing parameters

    # This function acts as an intermediate space to calculate new mass transfer coefficients
    # based on the linear relationship with alcohol. It also serves to calculate the generation/
    # depletion of each state variable

    v = types.SimpleNamespace()

    if type(t) == float:
        r = np.zeros([13, 1])
    else:
        r = np.zeros([13, t.size])

    # Define function for the pump overs switching on and off
    v.u = p.param_regression.u_max*(x.N/(p.param.K_N + x.N))*(x.S/(p.param.K_S2 + x.S))
    # Mass transfer coefficient for malic acid and relationship to ethanol concentration
    v.kla_mg_ma_0 = p.kla_mg_ma*(1 + p.param_regression.ema*x.E)

    # Term describing yeast death rate multiplied by biomass
    r[0, :] = x.X*(x.S/(p.param.K_S2 + x.S))
    # Term describing yeast specific growth rate multiplied by biomass
    r[1, :] = v.u*x.X

    # Term describing mass transfer of polymeric pigments from the must into the lees phase
    r[2, :] = p.param_regression.kla_ml_pp*(x.CPP - p.param_regression.PI_LSPP*(x.MPP_LS/x.M_LS))
    # Term reaction of anthocyanins to form polymeric pigments
    r[3, :] = p.param_regression.k_pp*x.CA

    # Term describing mass transfer of anthocyanins from the grapes into the must phase
    r[4, :] = p.kla_mg_a*(1 + p.param_regression.ea*x.E)*(p.param_regression.PI_GRA*x.CA_GR - x.CA)
    # Term describing mass transfer of anthocyanins from the must into the lees phase
    r[5, :] = p.param_regression.kla_ml_a*(x.CA - p.param_regression.PI_LSA*x.MA_LS/x.M_LS)

    # Term describing mass transfer of tannins from the grapes into the must phase
    r[6, :] = p.kla_mg_t*(1 + p.param_regression.et*x.E)*(p.param_regression.PI_GRT*x.CT_GR - x.CT)
    # Term describing mass transfer of tannins from the must into the lees phase
    r[7, :] = p.param_regression.kla_ml_t*(x.CT - p.param_regression.PI_LST*(x.MT_LS/x.M_LS))

    # Term describing mass transfer of tartaric acid from the grapes into the must phase
    r[8, :] = p.kla_mg_ta*(1 + p.param_regression.eta*x.E)*(p.param_regression.PI_GRTA*x.CTA_GR - x.CTA)
    # Term describing mass transfer of tartaric acid from the must into the lees phase
    r[9, :] = p.param_regression.kla_ml_ta*(x.CTA - p.param_regression.PI_LSTA*(x.MTA_LS / x.M_LS))

    # Term describing mass transfer of malic acid from the grapes into the must phase
    r[10, :] = p.kla_mg_ma*(1 + p.param_regression.ema*x.E)*(p.param_regression.PI_GRMA*x.CMA_GR - x.CMA)
    # Term describing mass transfer of malic acid from the must into the lees phase
    r[11, :] = p.param_regression.kla_ml_ma*(x.CMA - p.param_regression.PI_LSMA*(x.MMA_LS/x.M_LS))
    r[12, :] = p.param_regression.k_d*x.E*x.X

    # The next line of code converts the stoichiometric matrix into a dictionary in
    # order for it to be easily referenced when creating a new dictionary.

    SM = p.SM.__dict__

    # hold refers to a new dictionary, where the arrays contained in the stoichiometric
    # matrix are multiplied by the array containing the terms in array r
    hold = {p.svar_list[i]: SM[p.svar_list[i]].dot(r) for i in range(len(p.svar_list))}

    v.S = types.SimpleNamespace(**hold)

    return v


def ode(t, x_vec, p):
    # t - Time span of fermentation
    # u - A namespace of exogenous inputs - in this case the pump over regime
    # x - A name space of state variables used in conjunction with the stoichiometric matrix
    # p - A namespace containing parameters

    # This contains the system of ODEs for the fermentation/extraction model

    x = vector_to_namespace(x_vec, p.svar_list)
    v = int_var(t, x, p)

    ddt = types.SimpleNamespace()
    # The following equations describe the reactions and mass transfer occurring in the model

    ddt.X = v.S.X
    ddt.N = v.S.N
    ddt.S = v.S.S
    ddt.CO2 = v.S.CO2
    ddt.E = v.S.E
    ddt.GL = v.S.GL
    ddt.SC = v.S.SC
    ddt.CPP = v.S.CPP
    ddt.MPP_LS = v.S.MPP_LS
    ddt.CA = v.S.CA
    ddt.CA_GR = v.S.CA_GR
    ddt.MA_LS = v.S.MA_LS
    ddt.CT = v.S.CT
    ddt.CT_GR = v.S.CT_GR
    ddt.MT_LS = v.S.MT_LS
    ddt.CTA = v.S.CTA
    ddt.CTA_GR = v.S.CTA_GR
    ddt.MTA_LS = v.S.MTA_LS
    ddt.CMA = v.S.CMA
    ddt.CMA_GR = v.S.CMA_GR
    ddt.MMA_LS = v.S.MMA_LS
    ddt.M_LS = v.S.M_LS

    de = namespace_to_vector(ddt, p.svar_list)

    return de


def speed_up(p, x):
    # p:    A namespace of parameters
    # x:    A namespace of initial values for the model
    # The function splits up the model into different stages based on the occurrence of pump-overs, and declares
    # different maximum time steps in order to solve the model faster. The function will use the final solution values
    # of each stage as the initial values for the next stage and concatenate the solutions

    x_0 = namespace_to_vector(x, p.svar_list)
    my_array = np.zeros([len(x_0), 1])

    for i in range(1, 19):
        if i == 1:
            t = np.arange(0, i * 12 * 60)
            t2 = np.arange(i * 12 * 60, i * 12 * 60 + 6)

        elif i == 18:
            t = np.arange((i - 1) * 12 * 60 + 6, i * 12 * 60)
            t2 = np.arange(i * 12 * 60, i * 12 * 60 + 6)
            t3 = np.arange(i * 12 * 60 + 6, i * 12 * 60 + 60)

        else:
            t = np.arange((i - 1) * 12 * 60 + 6, i * 12 * 60)
            t2 = np.arange(i * 12 * 60, i * 12 * 60 + 6)

        p.kla_mg_a, p.kla_mg_t, p.kla_mg_ta, p.kla_mg_ma = p.param_regression.kla_mg_a_i, p.param_regression.kla_mg_t_i, p.param_regression.kla_mg_ta_i, p.param_regression.kla_mg_ma_i
        sol = solve_ivp(lambda t, x: ode(t, x, p), [t[0], t[-1]], x_0, method='LSODA', t_eval=t, max_step=60)
        x_0 = sol.y[:, -1]
        p.kla_mg_a, p.kla_mg_t, p.kla_mg_ta,  p.param.kla_mg_ma = p.param_regression.kla_max_a, p.param_regression.kla_max_t, p.param_regression.kla_max_ta, p.param_regression.kla_max_ma
        sol2 = solve_ivp(lambda t, x: ode(t, x, p), [t2[0], t2[-1]], x_0, method='LSODA', t_eval=t2, max_step=1)
        x_0 = sol2.y[:, -1]
        a = np.concatenate((sol.y, sol2.y), axis=1)

        if i == 18:
            p.kla_mg_a, p.kla_mg_t, p.kla_mg_ta, p.kla_mg_ma = p.param_regression.kla_mg_a_i, p.param_regression.kla_mg_t_i, p.param_regression.kla_mg_ta_i, p.param_regression.kla_mg_ma_i
            sol3 = solve_ivp(lambda t, x: ode(t, x, p), [t3[0], t3[-1]], x_0, method='LSODA', t_eval=t3, max_step=60)
            b = np.concatenate((a, sol3.y), axis=1)

        if i == 18:
            my_array = np.concatenate((my_array, b), axis=1)
        else:
            my_array = np.concatenate((my_array, a), axis=1)

    my_array = np.delete(my_array, 0, axis=1)
    return my_array


def simple_calc(t, x, p, time_array):
    # t - Time span of fermentation
    # x - A name space of state variables
    # p - A namespace containing parameters

    # This function acts as an intermediate space to calculate colour density, total phenolic index and pH
    # based on the values obtained during fermentation.

    K_W = 1.8e-16  # Ka value of water [unitless]
    K_SA = 6.21e-5  # Ka value of succinic acid [unitless]
    K_MA = 3.48e-4  # Ka value of malic acid [unitless]
    K_TA1 = 9.20e-4  # Ka value of tartaric acid with regard to first dissociation [unitless]
    K_TA2 = 4.31e-5  # Ka value of tartaric acid with regard to second dissociation [unitless]
    m_SC = (x.SC*20) / 118  # Total moles of succinic acid in 20L ferment
    m_MA = (x.CMA*20) / 134  # Total moles of malic acid in 20L ferment
    m_TA = (x.CTA*20) / 150  # Total moles of Tartaric acid in 20L ferment
    CI = 0.00000001  # Alkalinity of solution

    v = types.SimpleNamespace()
    # Must phase colour density
    v.cd = p.param_regression.a_CD * x.CA + p.param_regression.b_CD * x.CPP
    # Must phase total phenolic index
    v.tpi = p.param_regression.a_TPI1 * x.CA + p.param_regression.b_TPI1 * x.CPP + p.param_regression.c_TPI1 * x.CT
    # Cap phase total phenolic index
    v.tpi_gr = p.param_regression.a_TPI2 * x.CA_GR + p.param_regression.b_TPI2 * x.CT_GR

    pH = lambda h: K_W*h*(h + K_SA)*(h + K_MA)*(1 + K_TA1*h + K_TA1*K_TA2) \
                   + m_SC*(K_SA*h**2)*(h + K_SA)*(h + K_MA)*(1 + K_TA1*h + K_TA1*K_TA2) \
                   + m_MA*(K_MA*h**2)*(h + K_SA)*(h + K_MA)*(1 + K_TA1*h + K_TA1*K_TA2) \
                   + (m_TA*(K_TA1*h**2))*(h + K_SA)*(h + K_MA)*(1 + K_TA1*h + K_TA1*K_TA2)\
                   + (m_TA*(K_TA1*h) + 2*m_TA*K_TA1*K_TA2)*(h**2)*(h + K_SA)*(h + K_MA) \
                   - (CI*h**2 + h**3)*(h + K_SA)*(h + K_MA)*(1 + K_TA1*h + K_TA1 * K_TA2)

    pH_calc = []
    v.pH_must = []

    for i in range(0, len(time_array)):
        m_SC = (x.SC[time_array[i]]) / 118
        m_MA = (x.CMA[time_array[i]]) / 134
        m_TA = (x.CTA[time_array[i]]) / 150
        ans = fsolve(pH, 0.001)
        pH_calc.append(-np.log10(ans))

    for i in range(0, (9 * 60 * 24 + 60)):
        if i <= (1 * 24 * 60):
            v.pH_must.append(pH_calc[0])
        elif i < (1 * 24 * 60) and i <= (2 * 24 * 60):
            v.pH_must.append(pH_calc[1])
        elif i < (2 * 24 * 60) and i <= (3 * 24 * 60):
            v.pH_must.append(pH_calc[2])
        elif i < (3 * 24 * 60) and i <= (4 * 24 * 60):
            v.pH_must.append(pH_calc[3])
        elif i < (4 * 24 * 60) and i <= (5 * 24 * 60):
            v.pH_must.append(pH_calc[4])
        elif i < (5 * 24 * 60) and i <= (6 * 24 * 60):
            v.pH_must.append(pH_calc[5])
        elif i < (6 * 24 * 60) and i <= (7 * 24 * 60):
            v.pH_must.append(pH_calc[6])
        elif i < (8 * 24 * 60) and i <= (9 * 24 * 60):
            v.pH_must.append(pH_calc[7])
        else:
            v.pH_must.append(pH_calc[8])

    return v


def residuals(case, data, time_array, sd, p):
    # case:         A list of cases of form :['cultivar', 'temperature']
    # data:         A name space containing all experimental data (returned by experimental_data)
    # time_array:   An array of time points for the experimental data
    # sd:           A namespace of standard deviations (returned by experimental_data)
    # p:            A namespace of variables
    # The function returns a vector of experimental data points based on the case declared, as well as a vector of
    # standard deviations for the case declared

    exp_list = []
    stdev_list = []
    dict_sd = sd.__dict__
    dict_exp = {**data.Must.__dict__}

    case_types = types.SimpleNamespace(Specific=3, Temperature=9, Cultivar=9)

    if case[0] != 'n/a' and case[1] != 'n/a':
        det = case_types.Specific
    else:
        det = case_types.Temperature

    stdev_rec = np.zeros([len(p.interest * det), time_array.size])
    exp_values = np.zeros([len(p.interest * det), time_array.size])

    c = 0
    e = 0
    d = 0

    for i in range(0, (len(p.interest) * det)):
        if i % det == 0 and i != 0:
            c = c + 1
            if c == len(p.interest):
                c = len(p.interest)
        if c == 14:
            e = 0
        if c == 19:
            e = 0

        if i % 3 == 0 and i != 0:
            d = d + 1
            if i % det == 0:
                d = 0

        if case[0] != 'n/a' and case[1] != 'n/a':
            exp_list.append(
                case[0] + '_' + case[1] + '_' + 'Degree' + '_' + p.repeat[i % 3] + '_' + p.exp_components[c] + '_' +
                p.phase_name[e])

        if case[0] == 'n/a' and case[1] != 'n/a':
            exp_list.append( p.cultivars[d] + '_' + case[1] + '_' + 'Degree' + '_' + p.repeat[i % 3] + '_' +
                             p.exp_components[c] + '_' + p.phase_name[e])

        if case[0] != 'n/a' and case[1] == 'n/a':
            exp_list.append(case[0] + '_' + p.temperature[d] + '_' + 'Degree' + '_' + p.repeat[i % 3] + '_' +
                            p.exp_components[c] + '_' + p.phase_name[e])

        if case[0] != 'n/a' and case[1] != 'n/a':
            stdev_list.append(p.phase_name[e] + '_' + case[0] + '_' + case[1] + '_' + p.exp_components[c])

        if case[0] == 'n/a' and case[1] != 'n/a':
            stdev_list.append(p.phase_name[e] + '_' + case[1] + '_' + p.exp_components[c])

        if case[0] != 'n/a' and case[1] == 'n/a':
            stdev_list.append(p.phase_name[e] + '_' + case[0] + '_' + p.exp_components[c])

        exp_values[i, :] = dict_exp[exp_list[i]]
        stdev_rec[i] = dict_sd[stdev_list[i]]

    return exp_values.ravel(), stdev_rec.ravel(), det


def objective(p_vec, p, x, t, time_array, plist, det, exp_vec, stdev_vec, index):
    # p_vec:        A vector of initial guesses for the optimisation
    # p:            A namespace of parameters
    # x:            The initial values for the model
    # t:            A vector containing all the time points for the model
    # time_array:   A vector of time points for the experimental data
    # plist:        A list containing names of the vector of initial guesses
    # det:          The determinate returned by residuals, this declares the case being optimised
    # exp_vec:      A vector of experimental values (returned by residuals)
    # stdev_vec:    A vector of standard deviations (returned by residuals)
    # index:        A vector of integers, used to sort the predictions in the case of bootstrapping
    # This function calculates the solution to the model and finds the weighted difference between predicted values and
    # experimental values

    parameters = vector_to_namespace(p_vec, plist)
    print(parameters)
    p.param_regression = parameters

    a1 = speed_up(p, x)
    a11 = vector_to_namespace(a1, p.svar_list)
    a12 = simple_calc(t, a11, p, time_array)
    time_array[10] = 13019

    hold1 = {**a11.__dict__, **a12.__dict__}
    predictions = np.zeros([len(p.interest * det), time_array.size])
    c = 0
    pred_sample = []
    for i in range(0, (len(p.interest) * det)):
        if i % det == 0 and i != 0:
            c = c + 1
            if c == len(p.interest):
                c = len(p.interest)

        for j in range(0, len(time_array)):
            predictions[i, j] = hold1[p.interest[c]][time_array[j]]

    pred_hold = predictions.ravel()
    for i in range(0, len(index)):
        pred_sample.append(pred_hold[index[i]])
    pred = np.array(pred_sample)

    obj_hold = (exp_vec-pred)

    obj = np.zeros(len(obj_hold))
    for i in range(0, len(obj_hold)):
        obj[i] = obj_hold[i]/stdev_vec[i]

    return obj[~np.isnan(obj)]


def bill(p, det, time_array):
    # p:            A namespace of parameters
    # det:          The case (returned by residuals)
    # time_array:   A vector of time points for the experimental data
    # This is a function that uses sampling with replacement to generate a new vector of integers, used to sort the
    # predictions, experimental data, and standard deviations in the case of bootstrapping

    samples = np.arange(0, len(p.interest*det))
    bootstrap_index = random.choices(samples, k=len(p.interest*det*len(time_array)))

    return bootstrap_index


def davy_jones(p_vec, p, x, t, time_array, plist, det, exp_vec, stdev_vec, bounds, itt):
    # p_vec:        A vector of initial guesses for optimisation
    # p:            A namespace of parameters
    # x:            A namespaces of initial values for the model
    # t:            The timespan for the model
    # time_array:   A vector of time points for the experimental data
    # plist:        A list of names for the vector of initial guesses
    # det:          The case being optimised (returned by residuals)
    # exp_vec:      A vector of experimental values
    # stdev_vec:    A vector of standard deviations
    # bounds:       The bounds of the parameters to be optimised
    # itt:          The number of iterations for the bootstrapping
    # This function will perform multiple regressions on data that has een sample with replacement and return an array
    # of the parameters for each iteration

    x_array = np.zeros([itt, len(p.param_names_reg)])
    for i in range(0, itt):
        new_exp, new_sd = [], []
        btsp_index = bill(p, det, time_array)
        for j in range(0, len(btsp_index)):
            new_exp.append(exp_vec[btsp_index[i]]), new_sd.append([btsp_index[i]])

        a = least_squares(objective, p_vec, method='trf', jac='3-point',bounds=bounds, x_scale='jac',
                          tr_solver='lsmr', args=(p, x, t, time_array, plist, det, new_exp, new_sd,btsp_index),
                          verbose=2)
        x_array[i, :] = a.x

    return x_array


def subplots(ts, time_array, sol1, v1, sol2, v2, average, sd):
    # ts:           Time span of fermentation
    # time_array:   A vector of time points for the experimental data
    # sol1:         A name space of the solution for the unoptimised model
    # v1:           A name space containing calculated variables for the unoptimised model
    # sol2:         A name space of the solution for the optimised model
    # v2:           A name space containing calculated variables for the optimised model
    # average:      A name space containing the means of the experimental data
    # sd:           A name space containing the standard deviations of the experimental data
    # This function displays the graphs comparing optimised and unoptimised results

    fig, ((x, n, s), (co2, e, gl), (sc, mls, ph)) = plt.subplots(3, 3)
    fig.suptitle('Alcoholic Fermentation')
    x.set_title('Active Biomass')
    x.set(xlabel='time [min]', ylabel='Biomass [g/L]')
    x.plot(ts, sol1.X, 'g', label='Unoptimised')
    x.plot( ts, sol2.X, 'r',label='Optimised')
    x.errorbar(time_array, average.Must_20_Biomass,fmt='bo', yerr=sd.Must_20_Biomass, label= 'Experimental',
               elinewidth=1, capsize=2)
    x.legend()
    n.set_title('Nitrogen')
    n.set(xlabel='time [min]', ylabel='Nitrogen [g/L]')
    n.plot(ts, sol1.N, 'g', label='Unoptimised')
    n.plot(ts, sol2.N, 'r', label='Optimised')
    n.errorbar(time_array, average.Must_20_Nitrogen, fmt='bo', yerr=sd.Must_20_Nitrogen, label='Experimental',
               elinewidth=1, capsize=2)
    n.legend()
    s.set_title('Sugar')
    s.set(xlabel='time [min]', ylabel='Sugar [g/L]')
    s.plot(ts, sol1.S, 'g', label='Unoptimised')
    s.plot(ts, sol2.S, 'r', label='Optimised')
    s.errorbar(time_array, average.Must_20_Sugar, fmt='bo', yerr=sd.Must_20_Sugar, label='Experimental',
               elinewidth=1, capsize=2)
    s.legend()
    co2.set_title('Carbon dioxide')
    co2.set(xlabel='time [min]', ylabel='Carbon Dioxide [g/L]')
    co2.plot(ts, sol1.CO2, 'g', label='Unoptimised')
    co2.plot(ts, sol2.CO2, 'r', label='Optimised')
    co2.legend()
    e.set_title('Ethanol')
    e.set(xlabel='time [min]', ylabel='Ethanol [g/L]')
    e.plot(ts, sol1.E, 'g', label='Unoptimised')
    e.plot(ts, sol2.E, 'r', label='Optimised')
    e.errorbar(time_array, average.Must_20_Ethanol, fmt='bo', yerr=sd.Must_20_Ethanol, label='Experimental',
               elinewidth=1, capsize=2)
    e.legend()
    gl.set_title('Glycerol')
    gl.set(xlabel='time [min]', ylabel='Glycerol [g/L]')
    gl.plot(ts, sol1.GL, 'g', label='Unoptimised')
    gl.plot(ts, sol2.GL, 'r', label='Optimised')
    gl.legend()
    sc.set_title('Succinic Acid')
    sc.set(xlabel='time [min]', ylabel='Succinic Acid[g/L]')
    sc.plot(ts, sol1.SC, 'g', label='Unoptimised')
    sc.plot(ts, sol2.SC, 'r', label='Optimised')
    sc.errorbar(time_array, average.Must_20_Succinic_Acid, fmt='bo', yerr=sd.Must_20_Succinic_Acid, label='Experimental',
               elinewidth=1, capsize=2)
    sc.legend()
    mls.set_title('Lees Mass')
    mls.set(xlabel='time [min]', ylabel='Lees Mass [g]')
    mls.plot(ts, sol1.M_LS, 'g', label='Unoptimised')
    mls.plot(ts, sol2.M_LS, 'r', label='Optimised')
    mls.legend()
    ph.set_title('Ph')
    ph.set(xlabel='time [min]', ylabel='pH [-]')
    ph.plot(ts, v1.pH_must, 'g', label='Unoptimised')
    ph.plot(ts, v2.pH_must, 'r', label='Optimised')
    ph.errorbar(time_array, average.Must_20_pH, fmt='bo', yerr=sd.Must_20_pH, label='Experimental',
               elinewidth=1, capsize=2)
    ph.legend()

    fig2, ((cpp, ca, ct), (cta, cma, cd), (tpi, tpig, e2)) = plt.subplots(3, 3)
    fig.suptitle('Phenolics')
    cpp.set_title('Polymeric Pigments')
    cpp.set(xlabel='time [min]', ylabel='Polymeric Pigments [mg/L]')
    cpp.plot(ts, sol1.CPP, 'g', label='Unoptimised')
    cpp.plot(ts, sol2.CPP, 'r', label='Optimised')
    cpp.errorbar(time_array, average.Must_20_Polymeric_Pigments, fmt='bo', yerr=sd.Must_20_Polymeric_Pigments, label='Experimental',
               elinewidth=1, capsize=2)
    cpp.legend()
    ca.set_title('Anthocyanins')
    ca.set(xlabel='time [min]', ylabel='Anthocyanins [mg/L]')
    ca.plot(ts, sol1.CA, 'g', label='Unoptimised')
    ca.plot(ts, sol2.CA, 'r', label='Optimised')
    ca.errorbar(time_array, average.Must_20_Anthocyanins, fmt='bo', yerr=sd.Must_20_Anthocyanins,
                 label='Experimental',
                 elinewidth=1, capsize=2)
    ca.legend()
    ct.set_title('Tannins')
    ct.set(xlabel='time [min]', ylabel='Tannins [mg/L]')
    ct.plot(ts, sol1.CT, 'g', label='Unoptimised')
    ct.plot(ts, sol2.CT, 'r', label='Optimised')
    ct.errorbar(time_array, average.Must_20_Tannins, fmt='bo', yerr=sd.Must_20_Tannins,
                label='Experimental',
                elinewidth=1, capsize=2)
    ct.legend()
    cta.set_title('Tartaric Acid')
    cta.set(xlabel='time [min]', ylabel='Tartaric Acid [g/L]')
    cta.plot(ts, sol1.CTA, 'g', label='Unoptimised')
    cta.plot(ts, sol2.CTA, 'r', label='Optimised')
    cta.errorbar(time_array, average.Must_20_Tartaric_Acid, fmt='bo', yerr=sd.Must_20_Tartaric_Acid,
                label='Experimental',
                elinewidth=1, capsize=2)
    cta.legend()
    cma.set_title('Malic Acid')
    cma.set(xlabel='time [min]', ylabel='Malic Acid [g/L]')
    cma.plot(ts, sol1.CMA, 'g', label='Unoptimised')
    cma.plot(ts, sol2.CMA, 'r', label='Optimised')
    cma.errorbar(time_array, average.Must_20_Malic_Acid, fmt='bo', yerr=sd.Must_20_Malic_Acid,
                label='Experimental',
                elinewidth=1, capsize=2)
    cma.legend()
    cd.set_title('Colour Density')
    cd.set(xlabel='time [min]', ylabel='Colour Density [-]')
    cd.plot(ts, v1.cd, 'g', label='Unoptimised')
    cd.plot(ts, v2.cd, 'r', label='Optimised')
    cd.errorbar(time_array, average.Must_20_Colour_Density, fmt='bo', yerr=sd.Must_20_Colour_Density,
                label='Experimental',
                elinewidth=1, capsize=2)
    cd.legend()
    tpi.set_title('Total Phenolic Index')
    tpi.set(xlabel='time [min]', ylabel='Total Phenolic Index [-]')
    tpi.plot(ts, v1.tpi, 'g', label='Unoptimised')
    tpi.plot(ts, v2.tpi, 'r', label='Optimised')
    tpi.errorbar(time_array, average.Must_20_Total_Phenolic_Index, fmt='bo', yerr=sd.Must_20_Total_Phenolic_Index,
                label='Experimental',
                elinewidth=1, capsize=2)
    tpi.legend()
    tpig.set_title('Total Phenolic Index Cap')
    tpig.set(xlabel='time [min]', ylabel='Total Phenolic Index Cap [-]')
    tpig.plot(ts, v1.tpi_gr, 'g', label='Unoptimised')
    tpig.plot(ts, v2.tpi_gr, 'r', label='Optimised')
    tpig.errorbar(time_array, average.Cap_20_Total_Phenolic_Index, fmt='bo', yerr=sd.Cap_20_Total_Phenolic_Index,
                 label='Experimental',
                 elinewidth=1, capsize=2)
    tpig.legend()
    e2.set_title('Ethanol')
    e2.set(xlabel='time [min]', ylabel='Ethanol [g/L]')
    e2.plot(ts, sol1.E, 'g', label='Unoptimised')
    e2.plot(ts, sol2.E, 'r', label='Optimised')
    e2.errorbar(time_array, average.Must_20_Ethanol, fmt='bo', yerr=sd.Must_20_Ethanol, label='Experimental',
               elinewidth=1, capsize=2)
    e2.legend()

    fig3, ((cagr, ctgr), (ctagr, cmagr)) = plt.subplots(2, 2)
    fig.suptitle('Cap Phenolics')
    cagr.set_title('Anthocyanins')
    cagr.set(xlabel='time [min]', ylabel='Anthocyanins [mg/g]')
    cagr.plot(ts, sol1.CA_GR, 'g', label='Unoptimised')
    cagr.plot(ts, sol2.CA_GR, 'r', label='Optimised')
    cagr.errorbar(time_array, average.Cap_20_Anthocyanins, fmt='bo', yerr=sd.Cap_20_Anthocyanins, label='Experimental',
               elinewidth=1, capsize=2)
    cagr.legend()
    ctgr.set_title('Tannins')
    ctgr.set(xlabel='time [min]', ylabel='Tannins [mg/g]')
    ctgr.plot(ts, sol1.CT_GR, 'g', label='Unoptimised')
    ctgr.plot(ts, sol2.CT_GR, 'r', label='Optimised')
    ctgr.errorbar(time_array, average.Cap_20_Tannins, fmt='bo', yerr=sd.Cap_20_Tannins, label='Experimental',
                  elinewidth=1, capsize=2)
    ctgr.legend()
    ctagr.set_title('Tartaric Acid')
    ctagr.set(xlabel='time [min]', ylabel='Tartaric Acid [g/g]')
    ctagr.plot(ts, sol1.CTA_GR, 'g', label='Unoptimised')
    ctagr.plot(ts, sol2.CTA_GR, 'r', label='Optimised')
    ctagr.errorbar(time_array, average.Cap_20_Tartaric_Acid, fmt='bo', yerr=sd.Cap_20_Tartaric_Acid, label='Experimental',
                  elinewidth=1, capsize=2)
    ctagr.legend()
    cmagr.set_title('Malic Acid')
    cmagr.set(xlabel='time [min]', ylabel='Malic Acid [g/g]')
    cmagr.plot(ts, sol1.CMA_GR, 'g', label='Unoptimised')
    cmagr.plot(ts, sol2.CMA_GR, 'r', label='Optimised')
    cmagr.errorbar(time_array, average.Cap_20_Malic_Acid, fmt='bo', yerr=sd.Cap_20_Malic_Acid, label='Experimental',
                  elinewidth=1, capsize=2)
    cmagr.legend()

    fig4, ((mpp, ma, mt), (mta, mma, mls2)) = plt.subplots(2, 3)
    fig.suptitle('Lees Phenolics')
    mpp.set_title('Polymeric Pigments')
    mpp.set(xlabel='time [min]', ylabel='Polymeric Pigments[mg/L]')
    mpp.plot(ts, sol1.MPP_LS, 'g', label='Unoptimised')
    mpp.plot(ts, sol2.MPP_LS, 'r', label='Optimised')
    mpp.errorbar(time_array, average.Lees_20_Polymeric_Pigments, fmt='bo', yerr=sd.Lees_20_Polymeric_Pigments, label='Experimental',
               elinewidth=1, capsize=2)
    mpp.legend()
    ma.set_title('Anthocyanins')
    ma.set(xlabel='time [min]', ylabel='Anthocyanins[mg/L]')
    ma.plot(ts, sol1.MA_LS, 'g', label='Unoptimised')
    ma.plot(ts, sol2.MA_LS, 'r', label='Optimised')
    ma.errorbar(time_array, average.Lees_20_Anthocyanins, fmt='bo', yerr=sd.Lees_20_Anthocyanins,
                 label='Experimental',
                 elinewidth=1, capsize=2)
    ma.legend()
    mt.set_title('Tannins')
    mt.set(xlabel='time [min]', ylabel='Tannins[mg/L]')
    mt.plot(ts, sol1.MT_LS, 'g', label='Unoptimised')
    mt.plot(ts, sol2.MT_LS, 'r', label='Optimised')
    mt.errorbar(time_array, average.Lees_20_Tannins, fmt='bo', yerr=sd.Lees_20_Tannins,
                label='Experimental',
                elinewidth=1, capsize=2)
    mt.legend()
    mta.set_title('Tartaric Acid')
    mta.set(xlabel='time [min]', ylabel='Tartaric Acid[g/L]')
    mta.plot(ts, sol1.MTA_LS, 'g', label='Unoptimised')
    mta.plot(ts, sol2.MTA_LS, 'r', label='Optimised')
    mta.errorbar(time_array, average.Lees_20_Tartaric_Acid, fmt='bo', yerr=sd.Lees_20_Tartaric_Acid,
                label='Experimental',
                elinewidth=1, capsize=2)
    mta.legend()
    mma.set_title('Malic Acid')
    mma.set(xlabel='time [min]', ylabel='Malic Acid[mg/L]')
    mma.plot(ts, sol1.MMA_LS, 'g', label='Unoptimised')
    mma.plot(ts, sol2.MMA_LS, 'r', label='Optimised')
    mma.errorbar(time_array, average.Lees_20_Malic_Acid, fmt='bo', yerr=sd.Lees_20_Malic_Acid,
                label='Experimental',
                elinewidth=1, capsize=2)
    mma.legend()
    mls2.set_title('Lees Mass')
    mls2.set(xlabel='time [min]', ylabel='Lees Mass [g]')
    mls2.plot(ts, sol1.M_LS, 'g', label='Unoptimised')
    mls2.plot(ts, sol2.M_LS, 'r', label='Optimised')
    mls2.legend()
    plt.show()

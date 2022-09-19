# ======================================================================================================================
# Script Name: Main_Script.py
# File Name: Three_Phase_Final
# Date: 06-09-2022
# Creator: K.N. Lambrecht
# Purpose: This script contains the main parameters, and performs basic data handling
# ======================================================================================================================

import pandas as pd
from Function import *
import types
import numpy as np

# Define time and parameters
days = 9
minutes_fermentation = days*24*60 + 60.      # Calculates number of minutes that the fermentation runs for
ts = np.arange(0, minutes_fermentation)      # Define time space
V_MS = 20                                    # L, volume of tank
V_GR = 0.25*V_MS                             # Approximate volume of grapes
V_CL = 0.1*V_MS                              # Approximate volume of grapes

# # Implement Simple Name Space for the parameters of the system
# The units and the definitions are defined below

p = types.SimpleNamespace()     # Line included to time the optimisation

# This parameter space contains the variables which will not be changed by the regression
p.param = types.SimpleNamespace(K_S=9.843439862211241, K_S2=149.07302856445312, K_N=0.01,
                                K_E=1, Y_LX=0.00001, Y_NX=0.307, Y_SX=0.04,
                                Y_CX=0.29, Y_EX=0.28, Y_GX=0.14, Y_SCX=0.11,
                                B=0.0024352033212567746, Y_ES=3.441094938937915, Y_GS=0.003, Y_SCS=0.00001)

# This parameter space contains variables which will change with the regression
p.param_regression = types.SimpleNamespace(k_d=9.660648089497706e-06, u_max=0.0007970177076081221,  k_pp=5.500e-06,
                                           PI_LSPP=0.4, PI_GRA=0.76, PI_LSA=5, PI_GRT=0.36, PI_LST=1, PI_GRTA=0.36,
                                           PI_LSTA=1, PI_GRMA=2.14, PI_LSMA=1, a_CD=0.03, b_CD=0.02, a_TPI1=0.03,
                                           b_TPI1=0.02, c_TPI1=0.02, a_TPI2=0.01, b_TPI2=0.02, kla_ml_pp=0.01,
                                           kla_mg_a_i=0.0004, kla_ml_a=0.08, kla_mg_t_i=0.002, kla_ml_t=0.004,
                                           kla_mg_ta_i=0.0001, kla_ml_ta=0.0008, kla_mg_ma_i=0.00007, kla_ml_ma=0.0004,
                                           ea=0.002, et=0.0008460846885103096, eta=0.0001, ema=0.000002,
                                           kla_max_a=0.021010790130807933, kla_max_t=0.001,
                                           kla_max_ma=0.011010790130807933, kla_max_ta=0.011010790130807933)


# The lower and upper bounds for parameter selection are declared and transformed in vectors
lp = types.SimpleNamespace(k_d=9.5e-06, u_max=0.0005, k_pp=5.00e-06, PI_LSPP=0.0, PI_GRA=0, PI_LSA=1, PI_GRT=0,
                           PI_LST=0, PI_GRTA=0, PI_LSTA=0, PI_GRMA=0, PI_LSMA=0, a_CD=0.01, b_CD=0.01, a_TPI1=0.01,
                           b_TPI1=0.01, c_TPI1=0.01, a_TPI2=0.005, b_TPI2=0.01, kla_ml_pp=0.005, kla_mg_a_i=0.0002,
                           kla_ml_a=0.05, kla_mg_t_i=0.001, kla_ml_t=0.001, kla_mg_ta_i=0.00005, kla_ml_ta=0.0005,
                           kla_mg_ma_i=0.00005, kla_ml_ma=0.0001, ea=0.001, et=0.0005, eta=0.00005, ema=0.0000015,
                           kla_max_a=0.01, kla_max_t=0.0005, kla_max_ma=0.005, kla_max_ta=0.005)
#
up = types.SimpleNamespace(k_d=9.7e-06, u_max=0.001, k_pp=6.000e-06, PI_LSPP=0.5, PI_GRA=5, PI_LSA=5, PI_GRT=5,
                           PI_LST=1, PI_GRTA=5, PI_LSTA=1, PI_GRMA=5, PI_LSMA=1, a_CD=0.05, b_CD=0.05, a_TPI1=0.05,
                           b_TPI1=0.05, c_TPI1=0.05, a_TPI2=0.015, b_TPI2=0.05, kla_ml_pp=0.015, kla_mg_a_i=0.0005,
                           kla_ml_a=0.1, kla_mg_t_i=0.003, kla_ml_t=0.005, kla_mg_ta_i=0.0005, kla_ml_ta=0.001,
                           kla_mg_ma_i=0.0001, kla_ml_ma=0.0005, ea=0.005, et=0.001, eta=0.00015, ema=0.0000025,
                           kla_max_a=0.05, kla_max_t=0.002, kla_max_ma=0.015, kla_max_ta=0.015)


# The following is a stoichiometric matrix for the three-phase fermentation vessel

p.SM = types.SimpleNamespace(X=np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]),
                             N=np.array([0, -p.param.Y_NX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             S=np.array([-p.param.Y_ES*p.param.B, -p.param.Y_ES, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             CO2=np.array([p.param.B, p.param.Y_CX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             E=np.array([p.param.B, p.param.Y_EX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             GL=np.array([p.param.Y_GS, p.param.Y_GX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             SC=np.array([p.param.Y_SCS, p.param.Y_SCX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             CPP=np.array([0, 0, -1/V_MS, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             MPP_LS=np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                             CA=np.array([0, 0, 0, -1, 1/V_MS, -1/V_MS, 0, 0, 0, 0, 0, 0, 0]),
                             CA_GR=np.array([0, 0, 0, 0, -1/V_GR, 0, 0, 0, 0, 0, 0, 0, 0]),
                             MA_LS=np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]),
                             CT=np.array([0, 0, 0, 0, 0, 0, 1/V_MS, -1/V_MS, 0, 0, 0, 0, 0]),
                             CT_GR=np.array([0, 0, 0, 0, 0, 0, -1/V_GR, 0, 0, 0, 0, 0, 0]),
                             MT_LS=np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),
                             CTA=np.array([0, 0, 0, 0, 0, 0, 0, 0, 1/V_MS, -1/V_MS, 0, 0, 0]),
                             CTA_GR=np.array([0, 0, 0, 0, 0, 0, 0, 0, -1/V_GR, 0, 0, 0, 0]),
                             MTA_LS=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
                             CMA=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/V_MS, -1/V_MS, 0]),
                             CMA_GR=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/V_GR, 0, 0]),
                             MMA_LS=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),
                             M_LS=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, p.param.Y_LX*V_MS]))

# The following lists are used to index the variables of interest when the model is solved
p.svar_list = ['X', 'N', 'S', 'CO2', 'E', 'GL', 'SC', 'M_LS', 'CPP', 'MPP_LS', 'CA', 'CA_GR', 'MA_LS', 'CT', 'CT_GR',
               'MT_LS', 'CTA', 'CTA_GR', 'MTA_LS', 'CMA', 'CMA_GR', 'MMA_LS']

# The following lists contain the names of the parameters for indexing
p.param_names = ['K_S', 'K_S2', 'K_N', 'K_E', 'Y_LX', 'Y_NX', 'Y_SX', 'Y_CX', 'Y_EX', 'Y_GX', 'Y_SCX', 'B', 'Y_ES',
                 'Y_GS', 'Y_SCS']

p.param_names_reg = ['k_d', 'u_max', 'k_pp', 'PI_LSPP', 'PI_GRA', 'PI_LSA', 'PI_GRT', 'PI_LST', 'PI_GRTA', 'PI_LSTA',
                     'PI_GRMA', 'PI_LSMA', 'a_CD', 'b_CD', 'a_TPI1', 'b_TPI1', 'c_TPI1', 'a_TPI2', 'b_TPI2',
                     'kla_ml_pp', 'kla_mg_a_i', 'kla_ml_a', 'kla_mg_t_i', 'kla_ml_t', 'kla_mg_ta_i', 'kla_ml_ta',
                     'kla_mg_ma_i', 'kla_ml_ma', 'ea', 'et', 'eta', 'ema', 'kla_max_a', 'kla_max_t', 'kla_max_ma',
                     'kla_max_ta']

# The following lists are used to sort the data frames and build dictionaries surrounding the experimental data
p.scenario = ['All', 'Cultivar', 'Temperature', 'Specific']
p.cultivars = ['CS', 'SH', 'ME']
p.temperature = ['20', '25', '28']

p.cap_components = ['Anthocyanins', 'Tannins', 'Total_Phenolic_Index', 'Tartaric_Acid', 'Malic_Acid']
p.lees_components = ['Anthocyanins', 'Polymeric_Pigments', 'Tannins', 'Tartaric_Acid', 'Malic_Acid']
p.must_components = ['Anthocyanins', 'Polymeric_Pigments', 'Tannins', 'Total_Phenolic_Index', 'Tartaric_Acid',
                     'Malic_Acid', 'Colour_Density', 'Succinic_Acid', 'Biomass', 'Sugar', 'Ethanol', 'pH', 'CO2',
                     'Nitrogen']

p.repeat = ['Repeat_1', 'Repeat_2', 'Repeat_3']
p.phase_name = ['Must', 'Cap', 'Lees', 'Grape']
p.interest = ['X', 'N', 'S', 'CO2', 'E', 'SC', 'CPP', 'CA', 'CT', 'CTA', 'CMA', 'cd', 'tpi', 'pH_must']
p.exp_components = ['Biomass', 'Nitrogen', 'Sugar', 'CO2', 'Ethanol', 'Succinic_Acid', 'Polymeric_Pigments',
                    'Anthocyanins', 'Tannins', 'Tartaric_Acid', 'Malic_Acid', 'Colour_Density', 'Total_Phenolic_Index',
                    'pH']

# Import the experimental data
df_must = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\wine_data_N.csv')
df_lees = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\lees_data.csv')
df_cap = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\cap_data_normalised.csv')
df_grape = pd.read_csv('D:\\Stellenbosch\\2. PhD\\2022\\0. Experimental\\Data Sets\\grape_data.csv')

# Functions used to sort the data frames and obtain the averages and standard deviations
df = types.SimpleNamespace(Must=df_must, Lees=df_lees, Cap=df_cap, Grape=df_grape)
data, list_values = experimental_data(df, p.phase_name)
sd, average = standard_deviations(df, p.phase_name, p)

# This declares the time points for plotting experimental data
time_points = (df_must['Hour_Must']).to_numpy()
time_array = time_points*60
time_array[10] = 13019

# Case determines the type of optimisation. The first element is the cultivar, the second is the temperature of interest
# the options for the first element are 'CS', 'SH', 'ME', or 'n/a'; the options ofr the second element are '20', '25',
# '28' or 'n/a'
case1, case2, case3, case4, case5, case6, case7, case8, case9 = ['CS', '20'], ['CS', '25'], ['CS', '28'], ['SH', '20'], ['SH', '25'], ['SH', '28'], ['ME', '20'], ['ME', '25'], ['ME', '28']


# The following generates two vectors containing the experimental values and standard deviations for each point. Based
# on the case
exp_vector1, stdev_vector1, determinate1 = residuals(case1, data, time_array, sd, p)
exp_vector2, stdev_vector2, determinate2 = residuals(case2, data, time_array, sd, p)
exp_vector3, stdev_vector3, determinate3 = residuals(case3, data, time_array, sd, p)
exp_vector4, stdev_vector4, determinate4 = residuals(case4, data, time_array, sd, p)
exp_vector5, stdev_vector5, determinate5 = residuals(case5, data, time_array, sd, p)
exp_vector6, stdev_vector6, determinate6 = residuals(case6, data, time_array, sd, p)
exp_vector7, stdev_vector7, determinate7 = residuals(case7, data, time_array, sd, p)
exp_vector8, stdev_vector8, determinate8 = residuals(case8, data, time_array, sd, p)
exp_vector9, stdev_vector9, determinate9 = residuals(case9, data, time_array, sd, p)

# The following name space contains the initial values for the system of Ordinary Differential Equations
x1 = types.SimpleNamespace(X=0.25, N=average.Must_CS_20_Nitrogen[0], S=average.Must_CS_20_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_CS_20_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_CS_20_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_CS_20_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_CS_20_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_CS_20_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x2 = types.SimpleNamespace(X=0.25, N=average.Must_CS_25_Nitrogen[0], S=average.Must_CS_25_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_CS_25_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_CS_25_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_CS_25_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_CS_25_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_CS_25_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x3 = types.SimpleNamespace(X=0.25, N=average.Must_CS_28_Nitrogen[0], S=average.Must_CS_28_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_CS_28_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_CS_28_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_CS_28_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_CS_28_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_CS_28_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x4 = types.SimpleNamespace(X=0.25, N=average.Must_SH_20_Nitrogen[0], S=average.Must_SH_20_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_SH_20_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_SH_20_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_SH_20_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_SH_20_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_SH_20_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x5 = types.SimpleNamespace(X=0.25, N=average.Must_SH_25_Nitrogen[0], S=average.Must_SH_25_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_SH_25_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_SH_25_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_SH_25_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_SH_25_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_SH_25_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x6 = types.SimpleNamespace(X=0.25, N=average.Must_SH_28_Nitrogen[0], S=average.Must_SH_28_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_SH_28_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_SH_28_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_SH_28_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_SH_28_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_SH_28_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x7 = types.SimpleNamespace(X=0.25, N=average.Must_ME_20_Nitrogen[0], S=average.Must_ME_20_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_ME_20_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_ME_20_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_ME_20_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_ME_20_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_ME_20_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x8 = types.SimpleNamespace(X=0.25, N=average.Must_ME_25_Nitrogen[0], S=average.Must_ME_25_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_ME_25_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_ME_25_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_ME_25_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_ME_25_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_ME_25_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

x9 = types.SimpleNamespace(X=0.25, N=average.Must_ME_28_Nitrogen[0], S=average.Must_ME_28_Sugar[0], CO2=0, E=0, GL=0,
                           SC=0.015, M_LS=0.001, CPP=average.Must_ME_28_Polymeric_Pigments[0], MPP_LS=0.00,
                           CA=average.Must_ME_28_Anthocyanins[0], CA_GR=1009.6, MA_LS=0,
                           CT=average.Must_ME_28_Tannins[0], CT_GR=5384.49, MT_LS=0,
                           CTA=average.Must_ME_28_Tartaric_Acid[0], CTA_GR=3, MTA_LS=0,
                           CMA=average.Must_ME_28_Malic_Acid[0], CMA_GR=2.11, MMA_LS=0)

# Get lower and upper bounds into vectors and into a tuple
lower_bounds = namespace_to_vector(lp, p.param_names_reg)
upper_bounds = namespace_to_vector(up, p.param_names_reg)
bounds = (lower_bounds, upper_bounds)

# A new name space for parameters is generated for the optimisation
gs = p.param_regression
gl = p.param_names_reg
g = (namespace_to_vector(gs, gl))




# initial_p = [1, 1]
# a = np.linspace(0, 10, num=101)
# b1 = 0.5
# b2 = 5
# #
# normal = test_func(a, b1, b2)
# s = (test_func(a, b1, b2))*0.1
#
# numsam = 100
# noise = np.zeros((numsam, len(a)))
# for i in range(0, numsam):
#     rns = np.random.normal(size=len(a), scale=0.1)
#     noise[i, :] = normal + rns
#
# a_tile = np.tile(normal, (numsam, 1))
# sd_tile = np.tile(s, (numsam, 1))


# ======================================================================================================================
# Script Name: Bootstrap.py
# File Name: Three_Phase_Final
# Date: 13-09-2022
# Creator: K.N. Lambrecht
# Purpose: This script is used to perform the likelihood for the parameter sets
# ======================================================================================================================

import time
from mpire import WorkerPool
from matplotlib import pyplot as plt
import numpy as np
from Function import *
from Main_Script import *
from Bootstrap_variables import bootstrap_indices, num_iterations


if __name__ == '__main__':

    start = time.time()
    with WorkerPool(n_jobs=5) as pool:
        results = pool.map(regression, bootstrap_indices, iterable_len=num_iterations)

    stop = time.time()
    results_array = np.array(results)

    print(results_array, (stop-start))

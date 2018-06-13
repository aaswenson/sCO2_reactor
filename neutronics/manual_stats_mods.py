import numpy as np
import pandas as pd
import itertools
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.api as sms
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns
import plotting as plot

# apply linear or log shift to data
ops = {'lin' : lambda x: x,
       'log' : lambda x: np.log(x)
      }
data = pd.read_csv("depl_results.csv")
filter = ['keff > 1']
data = filter_data(filter, data)


def sm_lin_reg(predictors, target, form_predict='lin'):
    
    keys = [x[0] for x in predictors]
    X = data[keys]
    
    for var in predictors:
        X[var[0]] = ops[var[1]](X[var[0]])

    Y = ops[form_predict](data[target])
     
    model = sm.OLS(Y, X)
    results = model.fit()

    return results, model

predictors = [('core_r', 'log'), ('PD', 'log'), ('enrich', 'lin'), ('power', 'lin')]
res, model = sm_lin_reg(predictors, 'keff', 'lin')
print(res.summary())

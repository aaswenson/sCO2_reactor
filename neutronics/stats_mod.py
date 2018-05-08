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
data = filter_data([('keff', 'great', 1.0)], data)

def sk_lin_reg(predictors, target, train_frac, form_predict=('lin', 'lin')):

    X = ops[form_predict[0]](data[predictors])
    Y = ops[form_predict[1]](data[target])
    
    ntrain = int(len(X) * train_frac)
    X_train = X[:ntrain]
    X_test = X[ntrain:]

    Y_train = Y[:ntrain]
    Y_test = Y[ntrain:]

    regr = linear_model.LinearRegression()

    regr.fit(X_train, Y_train)

    y_pred = regr.predict(X_test)

    # The coefficients
    print('Coefficients: \n', regr.coef_)
    # The mean squared error
    print("Mean squared error: %.4f"
          % mean_squared_error(Y_test, y_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.4f' % r2_score(Y_test, y_pred))

    sse = np.sum(np.subtract(Y_test, y_pred)**2)

    print('The sum of squared errors is: %.4f' % sse)

    return r2_score(Y_test, y_pred), regr, sse


def sm_lin_reg(predictors, target, form_predict=('lin', 'lin')):

    X = ops[form_predict[0]](data[predictors])
    Y = ops[form_predict[1]](data[target])
    
    model = sm.OLS(Y, X)
    results = model.fit()

    return results, model

predictors = ['mass']
linlog = ('log', 'log')

res, model = sm_lin_reg(predictors, 'keff', linlog)
print(res.summary())
print("Now using sklearn linear regression.......")
sk_lin_reg(predictors, 'keff', 0.9, linlog)

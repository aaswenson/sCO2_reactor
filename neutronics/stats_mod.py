import numpy as np
import pandas as pd
import itertools
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.api as sms
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
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
    coef = np.ones(len(data))
    
    model = sm.OLS(Y, X)
    results = model.fit()

    return results, model

predictors = ['core_r', 'PD', 'enrich', 'power']#, 'cool_r']
linlog = ('lin', 'lin')

res, model = sm_lin_reg(predictors, 'keff', linlog)
print(res.summary())

## raw residuals vs. fitted
residsvfitted = plt.plot(res.predict(), res.resid,  'o')
l = plt.axhline(y = 0, color = 'grey', linestyle = 'dashed')
plt.xlabel('Fitted values')
plt.ylabel('Residuals')
plt.title('Residuals vs Fitted')
plt.show(residsvfitted)



#print("Now using sklearn linear regression.......")
#sk_lin_reg(predictors, 'keff', 0.9, linlog)

fig, ax = plt.subplots(figsize=(12,8))
fig = sm.graphics.plot_ccpr_grid(res, fig=fig)
#plt.show()

fig, ax = plt.subplots(figsize=(12,8))
fig = sm.graphics.influence_plot(res, ax=ax, criterion="cooks")
#plt.show()

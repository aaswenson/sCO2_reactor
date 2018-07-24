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
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns
import plotting as plot

# apply linear or log shift to data
ops = {'lin' : lambda x: x,
       'log' : lambda x: np.log(x),
       'sqr' : lambda x: np.power(x, 2)
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

def surf_plot(ind1, ind2, ind3, colorplot=None):
    """
    """

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X = ind1
    Y = ind2
    Z = ind3

    scaled_size = np.multiply(colorplot, 1000)

    # Plot the surface.
    cplt = colorplot
    p = ax.scatter(X,Y,Z, s=scaled_size, c=cplt,
            cmap=plt.cm.get_cmap('viridis',
        len(colorplot)))
    plt.colorbar(p, ax=ax, label='residuals')

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    plt.savefig('surface_plot.png', dpi=700)

    plt.show()

#predictors = [('core_r', 'lin'), ('PD', 'lin'), ('enrich', 'lin'), ('power',
#    'lin')]
predictors = [('enrich', 'log'), ('power', 'log'), ('core_r',
    'log'), ('PD', 'log')]
res, model = sm_lin_reg(predictors, 'keff', 'lin')
print(res.summary())

cmap_data = 'mass'
color_plot = data[cmap_data]
## raw residuals vs. fitted
#residsvfitted = plt.scatter(res.predict(), res.resid, c=color_plot,
#        cmap=plt.cm.get_cmap('plasma', len(set(color_plot))))
#plt.colorbar(label=cmap_data)
residsvfitted = plt.scatter(res.predict(), res.resid, s=3)
l = plt.axhline(y = 0, color = 'grey', linestyle = 'dashed')
plt.xlabel('Fitted values')
plt.ylabel('Residuals')
plt.title('Residuals vs Fitted')
plt.savefig('residuals_vs_fitted.png', dpi=700)
plt.show(residsvfitted)

#surf_plot(data['core_r'], data['PD'], data['enrich'], res.resid)

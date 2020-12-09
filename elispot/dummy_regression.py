# regression model of elispot assay
# example dataset from statsmodels module

# Generalized linear model 

# load modules
# import numpy as np
import statsmodels.api as sm

from scipy import stats
from matplotlib import pyplot as plt

from statsmodels import graphics
from statsmodels.graphics.api import abline_plot

# load the data and add a constant to the exogenous (independent) variables

# from config import config

print(sm.datasets.scotland.NOTE)

data = sm.datasets.scotland.load()

data.exog = sm.add_constant(data.exog)

# instatiate a model with default link function

gamma_model = sm.GLM(data.endog, data.exog, family=sm.families.Gamma())

# fit and summary

gamma_results = gamma_model.fit()

print(gamma_results.summary())

res = gamma_results

# quantities of interest

print('Total number of trials:',  data.endog[0].sum())
print('Parameters: ', res.params)
print('T-values: ', res.tvalues)

# plots

nobs = res.nobs
# y = data.endog[:,0]/data.endog.sum()
y = data.endog[:]/data.endog.sum()
yhat = res.mu

fig, ax = plt.subplots()
ax.scatter(yhat, y)
line_fit = sm.OLS(y, sm.add_constant(yhat, prepend=True)).fit()
abline_plot(model_results=line_fit, ax=ax)


ax.set_title('Model Fit Plot')
ax.set_ylabel('Observed values')
ax.set_xlabel('Fitted values')

# plot yhat vs y
fig, ax = plt.subplots()
ax.scatter(yhat, y)
line_fit = sm.OLS(y, sm.add_constant(yhat, prepend=True)).fit()
abline_plot(model_results=line_fit, ax=ax)


ax.set_title('Model Fit Plot')
ax.set_ylabel('Observed values')
ax.set_xlabel('Fitted values')

# plot yhat vs. Pearson residuals

fig, ax = plt.subplots()

ax.scatter(yhat, res.resid_pearson)
ax.hlines(0, 0, 1)
ax.set_xlim(0, 1)
ax.set_title('Residual Dependence Plot')
ax.set_ylabel('Pearson Residuals')
ax.set_xlabel('Fitted values')

# histogram of standardized deviance residuals

fig, ax = plt.subplots()

resid = res.resid_deviance.copy()
resid_std = stats.zscore(resid)
ax.hist(resid_std, bins=25)
ax.set_title('Histogram of standardized deviance residuals')

# QQ Plot of deviance residuals

graphics.gofplots.qqplot(resid, line='r')

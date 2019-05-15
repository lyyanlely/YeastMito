# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 09:52:51 2017

@author: kuns
"""

from __future__ import division
import numpy as np


def mean(data):
    '''
    Computes mean of an itemfreq array.
    If cdat = itemfreq(dat), then mean(cdat) = np.mean(dat)
    '''
    norm = np.sum(data[:,1])
    data_vals = data[:,0]
    data_nums = data[:,1]
    return np.sum(data_vals * data_nums) / norm


def var(data, ddof=0):
    '''
    Computes variance of an itemfreq array.
    If cdat = itemfreq(dat), then var(cdat) = np.var(dat)
    ddof: delta degrees of freedom (same as np)
    '''
    norm = np.sum(data[:,1]) - ddof
    data_vals = data[:,0]
    data_nums = data[:,1]
    m = mean(data)
    return np.sum((data_vals - m)**2 * data_nums) / norm


def std(data, ddof=0):
    '''
    Computes standard deviation of an itemfreq array.
    If cdat = itemfreq(dat), then std(cdat) = np.std(dat)
    ddof: delta degrees of freedom (same as np)
    '''
    return np.sqrt(var(data, ddof=ddof))


def skew(data):
    '''
    Computes the skewness of an itemfreq array.
    Note that this is the Pearson's moment coefficient of skewness which
    is the third central moment.
    If cdat = itemfreq(dat), then skew(cdat) = stats.skew(dat)
    '''
    norm = np.sum(data[:,1])
    data_vals = data[:,0]
    data_nums = data[:,1]
    m = mean(data)
    m3 = np.sum((data_vals - m) ** 3 * data_nums) / norm
    m2 = np.sum((data_vals - m) ** 2 * data_nums) / norm
    return m3 / (m2 ** 1.5)


def kurtosis(data):
    '''
    Computes the curtosis of an itemfreq array.
    Note that this is the fourth central moment.
    If cdat = itemfreq(dat), then kurtosis(cdat) = stats.kurtosis(dat)
    '''
    norm = np.sum(data[:,1])
    data_vals = data[:,0]
    data_nums = data[:,1]
    m = mean(data)
    m4 = np.sum((data_vals - m) ** 4 * data_nums) / norm
    m2 = np.sum((data_vals - m) ** 2 * data_nums) / norm
    return m4 / (m2 ** 2)


def median(data):
    '''
    Computes the median of an itemfreq array.
    If cdat = itemfreq(dat), then median(cdat) = stats.median(dat)
    '''
    return quantile(data, 0.5)


def quantile(data, q):
    '''
    Computes the q-quantile of an itemfreq array.
    I.e. returns the x for which
    P(X < x) = q
    '''
    x, c = cdf(data)
    low = c < q
    try:
        return x[np.nonzero(low)[0][-1]]
    except IndexError:
        # more than 1/q of the data has the lowest value
        return x[0]
            

def cdf(data):
    '''
    Computes the empirical cdf from an itemfreq array.
    (plt.hist(data[:,0], weights=data[:,1]) would give
    the pdf of the cdf computed by this function)
    then plt.plot(x, cdf) gives you the cdf.
    '''
    cdf = np.cumsum(data[:,1]) / np.sum(data[:,1])
    return data[:,0], cdf
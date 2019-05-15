# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 12:18:15 2017

@author: kuns
"""

from __future__ import division
import numpy as np
from collections import Counter


def mean(data, invert=False):
    '''
    Computes mean of a counter.
    If cdat = Counter(dat), then mean(cdat) = np.mean(dat)
    '''
    if invert:
        data = Counter({1/k: v for k, v in data.items()})
    norm = np.sum(data.values())
    data_vals = np.array(data.keys())
    data_nums = np.array(data.values())
    return np.sum(data_vals * data_nums) / norm


def var(data, ddof=0, invert=False):
    '''
    Computes variance of a counter.
    If cdat = Counter(dat), then var(cdat) = np.var(dat)
    ddof: delta degrees of freedom (same as np)
    '''
    if invert:
        data = Counter({1/k: v for k, v in data.items()})
    norm = np.sum(data.values()) - ddof
    data_vals = np.array(data.keys())
    data_nums = np.array(data.values())
    m = mean(data)
    return np.sum((data_vals - m) ** 2 * data_nums) / norm


def std(data, ddof=0, invert=False):
    '''
    Computes standard deviation of a counter.
    If cdat = Counter(dat), then std(cdat) = np.std(dat)
    ddof: delta degrees of freedom (same as np)
    '''
    return np.sqrt(var(data, ddof=ddof, invert=invert))


def skew(data):
    '''
    Computes the skewness of a counter.
    Note that this is the Pearson's moment coefficient of skewness which
    is the third central moment.
    If cdat = counter(dat), then skew(cdat) = stats.skew(dat)
    '''
    norm = np.sum(data.values())
    data_vals = np.array(data.keys())
    data_nums = np.array(data.values())
    m = mean(data)
    m3 = np.sum((data_vals - m) ** 3 * data_nums) / norm
    m2 = np.sum((data_vals - m) ** 2 * data_nums) / norm
    return m3 / (m2 ** 1.5)


def kurtosis(data):
    '''
    Computes the curtosis of a counter.
    Note that this is the fourth central moment.
    If cdat = counter(dat), then kurtosis(cdat) = stats.kurtosis(dat)
    '''
    norm = np.sum(data.values())
    data_vals = np.array(data.keys())
    data_nums = np.array(data.values())
    m = mean(data)
    m4 = np.sum((data_vals - m) ** 4 * data_nums) / norm
    m2 = np.sum((data_vals - m) ** 2 * data_nums) / norm
    return m4 / (m2 ** 2)


def median(data, invert=False):
    '''
    Computes the median of a counter.
    If cdat = counter(dat), then median(cdat) = stats.median(dat)
    '''
    return quantile(data, 0.5, invert=invert)


def quantile(data, q, invert=False):
    '''
    Computes the q-quantile of a counter
    I.e. returns the x for which
    P(X < x) = q
    '''
    x, c = cdf(data, invert=invert)
    low = c < q
    try:
        return x[np.nonzero(low)[0][-1]]
    except IndexError:
        # More than 1/q of the data has the lowest value
        return x[0]


def cdf(data, norm=True, invert=False):
    '''
    Computes the empirical cdf from a counter with values the weights of keys.
    (plt.hist(data.keys(), weights=data.values()) would give
    the pdf of the cdf computed by this function)
    then plt.plot(x, cdf) gives you the cdf.
    '''
    if invert:
        data = Counter({1/k: v for k, v in data.items()})
    sorted_data = sorted(data.items())
    x = np.array([y[0] for y in sorted_data])
    counts = np.array([y[1] for y in sorted_data])
    if norm:
        cdf = np.cumsum(counts) / np.sum(counts)
    else:
        cdf = np.cumsum(counts)
    return x, cdf
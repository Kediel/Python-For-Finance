#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:33:31 2017

@author: KDL
"""
#
#Valuation of European call options in Black-Scholes-Merton model
#including Vega function and implied votality estimation
#bsm_functions.py
#
from math import log, sqrt, exp
from scipy import stats
#Analytical BSM Formula

def bsm_call_value(S0, K, T, r, sigma):
    "Valuation of European call option in BSM model. Parameter definitions included on page 52 of Python in Finance"
    
    S0 = float(S0)
    d1 = (log(S0 / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * sqrt(T))
    d2 = (log(S0 / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * sqrt(T))
    value = (S0 * stats.norm.cdf(d1, 0.0, 1.0) - K * exp(-r * T) * stats.norm.cdf(d2, 0.0, 1.0))
    #stats.norm.cdf --> cumulative distribution function for normal distribution
    return value

#Vega Function

def bsm_vega(S0, K, T, r, sigma):
    "Vega of European option in BSM model."
    
    S0 = float(S0)
    d1 = (log(S0 / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * sqrt(T))
    vega = S0 * stats.norm.cdf(d1, 0.0, 1.0) * sqrt(T)
    return vega 
    
#Implied votality function

def bsm_call_imp_vol(S0, K, T, r, C0, sigma_est, it=100):
    "Implies volatility of European call option in BSM model"

    for i in range(it):
        sigma_est -= ((bsm_call_value(S0, K, T, r, sigma_est) - C0) / bsm_vega(S0, K, T, r, sigma_est))
    return sigma_est
    
    

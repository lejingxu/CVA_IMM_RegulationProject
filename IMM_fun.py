import numpy as np
from scipy.stats import norm


def CalcEEPE(EE, freq):
    '''
    Calculate Effective Expective Positive Exposure from Expectve Exposure profile
    Args:
        freq: simulation frequency per year. Ex. 12 for monthly grid, 52 for weekly, 252 for daily
    '''
    EE_runningMax = np.maximum.accumulate(EE[:freq])
    return np.mean(EE_runningMax)


def CalcEAD(Basel, EEPE, CVA, alpha, freq):
    '''
    Calculate 1 year Exposure at Default using Expective Exposure profile, according to Basel II or III
    Args:
        freq: simulation frequency per year. Ex. 12 for monthly grid, 52 for weekly, 252 for daily
    '''
    #EEPE = CalcEEPE(EE, freq)
    if Basel == '2':
        EAD = alpha * EEPE
    elif Basel == '3':
        EAD = 1.06 * np.maximum(alpha * EEPE - CVA, 0)
    else:
        EAD = 0
        print('wrong Basel accord number')
    return EAD


def EffectiveMaturity_old(EE, freq):
    '''
    Calculate effective maturity
    Args:
        freq: simulation frequency per year. Ex. 12 for monthly grid, 52 for weekly, 252 for daily
    '''
    EE_runningMax = np.maximum.accumulate(EE[:freq])
    M = np.maximum(1, np.minimum(5, (np.sum(EE_runningMax + np.sum(EE[freq:]))) / np.sum(EE_runningMax)))
    
    return M
    
def EffectiveMaturity(EE, freq,P_OISs):
    '''
    Calculate effective maturity
    Args:
        freq: simulation frequency per year. Ex. 12 for monthly grid, 52 for weekly, 252 for daily
    '''
    EE_runningMax = np.maximum.accumulate(EE[:freq])
    #M = np.maximum(1, np.minimum(5, (np.sum(EE_runningMax + np.sum(EE[freq:]))) / np.sum(EE_runningMax)))
    mean_EEE = np.average(np.multiply(EE_runningMax,P_OISs[0][:freq]))
    mean_EE = np.average(np.multiply(EE[freq:],P_OISs[0][freq:]))
    M = np.maximum(1, np.minimum(5,(mean_EEE+mean_EE)/mean_EEE))
    return M


def RegulatoryCapital_old(EAD, M, LGD, PD):
    '''
    Calculate regulatory capital
    Args:
    '''
    q = 0.001
    b = (0.11852-0.05478*np.log(PD))**2
    k = (1 + (M - 2.5) * b) / (1 - 1.15 * b)
    rho = 0.24 - 0.12 * (1-np.exp(-50*PD))
    RW = LGD * norm.cdf((norm.ppf(PD) - (np.sqrt(rho) * norm.ppf(q))) / np.sqrt(1. - rho))*k
    RC = EAD * RW
    return RC
    
def RegulatoryCapital(EAD, M, LGD, PD):
    '''
    Calculate regulatory capital
    Args:
    '''
    q = 0.001
    b = (0.11852-0.05478*np.log(PD))**2
    k = (1. + (M - 2.5) * b) / (1. - 1.5 * b)
    rho = 0.24 - 0.12 * (1-np.exp(-50.*PD))
    RW = LGD * norm.cdf((norm.ppf(PD) - (np.sqrt(rho) * norm.ppf(q))) / np.sqrt(1. - rho))*k
    RC = EAD * RW
    return RC
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 03 10:25:19 2017

@author: Lejing
"""

import numpy as np

def calculatePVEE(lbda,P_OIS,X,V_df,switch_collateral,switch_downProv,collateral,D):
    '''
    Args:
        ** Note that all simulated variables are of size num_simulation, each 
        ** element is the variable for one simulation
        lbda:              default intensity
        P_OIS:             discount bound price using OIS
        X:                 survival probability
        V_df:              default-free value of the portfolio
        switch_collateral: True or False for collateral agreement
        switch_downProv:   True or False for downgrade provision
        collateral:        collateral threshold
        D:                 default intensity threshold
        
    Return:
        
    '''
    num_simulation = len(V_df)
    #print "numsim",num_simulation
    PVEE = []
    EE = []
    
    exposures = [calculateExposure(V_df[k],switch_collateral,switch_downProv,collateral,D,lbda[k]) for k in range(num_simulation)]
    #print "exposures",exposures 
    for i in range(len(V_df[0])):
        numer = np.zeros(num_simulation)
        denom = np.zeros(num_simulation)
        numer_ee = np.zeros(num_simulation)
        for j in range(num_simulation):
            denom[j] = X[j][0][i]*lbda[j][i][0]
            numer_ee[j] = denom[j]*exposures[j][i]
            numer[j] = numer_ee[j]*P_OIS[j][0][i]
            
    
        numerator = np.average(numer)
        denominator = np.average(denom)
        numerator_ee = np.average(numer_ee)
        pvee = numerator/denominator
        ee = numerator_ee/denominator
        
        PVEE.append(pvee)
        EE.append(ee)
    
    return PVEE,EE

def calculateExposure(V_df,switch_collateral,switch_downProv,collateral,D,lbda):
    '''
    Args:
       
        
    Return:
        
    '''
    if not (switch_collateral or switch_downProv):
        return np.maximum(V_df,0)
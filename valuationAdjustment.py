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
        
    if switch_downProv:
        # extract instantaneous default intensity lbda(t_i,t_(i+1))
        lbda_inst = [lbda[i][0] for i in range(1,len(lbda),1)] 
        terminationTime = np.where(lbda_inst>D)[0]
        V_downProv = np.asarray(V_df)
        if terminationTime!= []:
            V_downProv[terminationTime[0]:] = 0
        EE = np.maximum(V_downProv,0)
        
        if switch_collateral:
            return np.maximum(np.minimum(EE,collateral),0)
        else:
            return EE
    else:
        return np.maximum(np.minimum(V_df,collateral),0)

        
def calculateUniCVA(EE,P_OIS,X,lbda,rr):
    '''
    Args:
        ** Note that all simulated variables are of size num_simulation, each 
        ** element is the variable for one simulation
        lbda:              default intensity
        P_OIS:             discount bound price using OIS
        X:                 survival probability
        rr:                recovery rate
        EE:                expected exposure
    '''
    CVA = (1.-rr)*np.sum(np.multiply(np.multiply(EE,P_OIS[0][0]),np.multiply(X[0][0],lbda[0][0])))
    return CVA
    
def calculateUniDVA(EE,P_OIS,X,lbda,rr):
    return calculateUniCVA(EE,P_OIS,X,lbda,rr)
    
def calculateNetUniCVA(uniCVA,uniDVA):
    return uniDVA - uniCVA
    
def calculateBiCVA(EE,P_OIS,X,lbda,rr,X_self):
    biCVA = (1.-rr)*np.sum(np.multiply(np.multiply(np.multiply(EE,X_self[0][0]),P_OIS[0][0]),np.multiply(X[0][0],lbda[0][0])))
    return biCVA
    
def calculateBiDVA(EE,P_OIS,X,lbda,rr,X_self):
    return calculateBiCVA(EE,P_OIS,X,lbda,rr,X_self)
    
def calculateNetBiCVA(biCVA,biDVA):
    return biCVA-biDVA
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 02 13:47:25 2017

@author: Lejing
"""
import numpy as np

def simulateXt_1(ts,sigma,kappa,wt,ind):
    '''
    One-factor Gaussian Model
    Args:
        ind: index of which row of wt to use, 2 for lbdaB, 3 lbdaC
        
    '''
    yt = sigma*sigma/2./kappa*(1.-np.exp(-2.*kappa*ts))
    xt = []
    x = 0
    dt = ts[1]-ts[0]
    for i in range(len(ts)):
        x = x + (yt[i] - kappa*x)*dt+sigma*wt[ind][i]
        xt.append(x)
        
    return xt,yt
    


def simulateXt(rho_x, sigma1, sigma2, kappa1, kappa2, maturity,ts,wt):
    '''
    Two-factor Gaussian model
    Args:
        
    
    '''
    sigma22 = sigma2
    sigma21 = sigma1*rho_x
    sigma11 = np.sqrt(sigma1**2 *(1 - rho_x*rho_x))
    
    dt = ts[1]-ts[0]
    x = np.asarray([0,0])
    xt = []
    yt = []
    
    for i in range(int(maturity/dt)):
        B = np.diag([np.exp(-kappa1*ts[i]),np.exp(-kappa2*ts[i])])
        
        integral_aa11 = (sigma11**2+sigma21**2)/2./kappa1*(np.exp(2*kappa1*ts[i])-1.)
        integral_aa12 = sigma21*sigma22/(kappa1+kappa2)*(np.exp((kappa1+kappa2)*ts[i])-1.)
        integral_aa22 = sigma22**2 /2./kappa2*(np.exp(2*kappa2*ts[i]-1.))
        integral_aa = np.asarray([[integral_aa11,integral_aa12],[integral_aa12,integral_aa22]])
        
        y = B.dot(integral_aa).dot(B)
        wi = np.array([wt[2][i],wt[3][i]])
        
        x = x + (y.dot(np.ones(2))-np.diag([kappa1,kappa2]).dot(x))*dt+np.diag([sigma1,sigma2]).dot(wi)
        xt.append(x)
        yt.append(y)
        
    return xt,yt
    
    
def simulateOIS(rho_x, sigma1, sigma2, kappa1, kappa2, sim_freq, maturity, f0_OIS, spread,ts,Tis,wt):
    
    P_OIS = []
    P_LIBOR = []
    
    P_OIS.append([np.exp(-f0_OIS*t) for t in ts]) # calculate P_OIS(0,t) for all t
    P_LIBOR.append([P_OIS[0][i]*np.exp(-spread*ts[i]) for i in range(ts.shape[0])])# calculate P_LIBOR(0,t) for all t
    xt,yt = simulateXt(rho_x, sigma1, sigma2, kappa1, kappa2, maturity,ts,wt) #120 xt's and yt's
    
    dt = ts[1]-ts[0]
    dT = Tis[1]-Tis[0]

    for i in range(ts.shape[0]):
        P_OIS_t = []
        P_OIS_0_t = P_OIS[0][i]
        P_LIBOR_t = []
        for j in np.where(Tis>=ts[i])[0]: # index for all coupon paying dates after t
            G = np.asarray([(1.-kappa1*(Tis[j]-ts[i]))/kappa1,(1.-kappa2*(Tis[j]-ts[i]))/kappa2])
            P_OIS_0_Ti = P_OIS[0][int((j+1)*dT/dt-1)]
            P_OIS_t_Ti = P_OIS_0_Ti/P_OIS_0_t*np.exp(-G.T.dot(xt[i])-0.5*G.T.dot(yt[i]).dot(G))
            P_OIS_t.append(P_OIS_t_Ti)
            P_LIBOR_t.append(P_OIS_t_Ti*np.exp(-spread*(Tis[j]-ts[i])))
            
        P_OIS.append(P_OIS_t)
        P_LIBOR.append(P_LIBOR_t)
        
        
    #print P_LIBOR
        
    return P_OIS,P_LIBOR
    
    
    
def simulateOIS_new(rho_x, sigma1, sigma2, kappa1, kappa2, sim_freq, maturity, f0_OIS, spread,ts,wt):
    '''
    *** DO NOT USE THIS FUNCTION. NOT FIXED. USE simulateOIS.
    '''
    
    P_OIS = []
    P_LIBOR = []
    
    P_OIS.append([np.exp(-f0_OIS*t) for t in ts]) # calculate P_OIS(0,t) for all t
    P_LIBOR.append([P_OIS[0][i]*np.exp(-spread*ts[i]) for i in range(ts.shape[0])])# calculate P_LIBOR(0,t) for all t
    xt,yt = simulateXt(rho_x, sigma1, sigma2, kappa1, kappa2, maturity,ts,wt) #120 xt's and yt's
    

    for i in range(ts.shape[0]):
        P_OIS_t = []
        P_LIBOR_t = []
        P_OIS_0_t = P_OIS[0][i]
        
        for j in np.where(ts>=ts[i])[0]: 
            G = np.asarray([(1.-kappa1*(ts[j]-ts[i]))/kappa1,(1.-kappa2*(ts[j]-ts[i]))/kappa2])
            P_OIS_0_tj = P_OIS[0][j]
            P_OIS_t_tj = P_OIS_0_tj/P_OIS_0_t*np.exp(-G.T.dot(xt[i])-0.5*G.T.dot(yt[i]).dot(G))
            P_OIS_t.append(P_OIS_t_tj)
            P_LIBOR_t.append(P_OIS_t_tj*np.exp(-spread*(ts[j]-ts[i])))
            
        P_OIS.append(P_OIS_t)
        P_LIBOR.append(P_LIBOR_t)
        
    return P_OIS,P_LIBOR    
    
def simulateSurvivalProb(lbda0_B,lbda0_C,ts,sigmaB,kappaB,sigmaC,kappaC,wt):
    XB = []
    XC = []
    lbdaB = []
    lbdaC = []

    XB.append([np.exp(-lbda0_B*t) for t in ts])
    XC.append([np.exp(-lbda0_C*t) for t in ts])
    lbdaB.append([lbda0_B]*len(ts))
    lbdaC.append([lbda0_C]*len(ts))
    
    xt_B,yt_B = simulateXt_1(ts,sigmaB,kappaB,wt,0)
    xt_C,yt_C = simulateXt_1(ts,sigmaC,kappaC,wt,1)
    
    for i in range(ts.shape[0]):
        XB_ti = []
        XC_ti = []
        lbdaB_ti = []
        lbdaC_ti = []
        for j in np.where(ts>=ts[i])[0]:
            G_B = (1.-np.exp(-kappaB*(ts[j]-ts[i]))/kappaB)
            XB_ti_tj = XB[0][j]/XB[0][i]*np.exp(-xt_B[i]*G_B-0.5*yt_B[i]*G_B*G_B)
            lbdaB_ti_tj = lbdaB[0][j]+np.exp(-kappaB*(ts[j]-ts[i]))*(xt_B[i]+yt_B[i]*G_B)
            XB_ti.append(XB_ti_tj)
            lbdaB_ti.append(lbdaB_ti_tj)
            
            G_C = (1.-np.exp(-kappaC*(ts[j]-ts[i]))/kappaC)
            XC_ti_tj = XC[0][j]/XC[0][i]*np.exp(-xt_C[i]*G_C-0.5*yt_C[i]*G_C*G_C)
            lbdaC_ti_tj = lbdaC[0][j]+np.exp(-kappaC*(ts[j]-ts[i]))*(xt_C[i]+yt_C[i]*G_C)
            XC_ti.append(XC_ti_tj)
            lbdaC_ti.append(lbdaC_ti_tj)
            
        XB.append(XB_ti)
        XC.append(XC_ti)
        lbdaB.append(lbdaB_ti)
        lbdaC.append(lbdaC_ti)
        
    return XB,XC,lbdaB,lbdaC
            
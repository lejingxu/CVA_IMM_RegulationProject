# -*- coding: utf-8 -*-
"""
Created on Sat Dec 02 10:25:15 2017

@author: Lejing
"""
import numpy as np


class Swap(object) :
    def __init__(self, maturity, coupon, freq, notional) :
        '''
        Args:
            maturity: in years
            coupon:   fixed coupon, in decimal
            freq:     frequency of payments, 2 for semi-annual
        
        '''
        self.maturity = maturity
        self.coupon = coupon
        self.freq = freq
        self.notional = notional
        
    def __str__(self) :
        return 'Swap: maturity %g, coupon %g, freq %g, notional %g'% \
        (self.maturity, self.coupon, self.freq, self.notional)
    
def priceSwap(swap, payerOrReceiver, P_OIS, P_LIBOR, i, ts, Tis, sim_freq) :
    '''
    Args:
        swap:     a swap object
        payerOrReceiver: string type, 'payer' or 'receiver'
        P_OIS:    array, discount bond price using OIS
        P_Libor:  array, discount bond price using LIBOR
        i:        time step now, takes positive integer
        ts:       simulation time steps
        Tis:      payment time steps
        
    return:
        price:    value of swap at time ts[i]
    '''
    #c = swap.coupon * swap.notional/swap.freq # fixed coupon
    c = swap.coupon/swap.freq
    dT = Tis[1]-Tis[0]
    startT = np.where(Tis>ts[i])[0]
    floating = []
    
    if len(startT) != 0:
        # add the first stub payment
        first_stub = startT[0]
        if first_stub != 0:
            floating.append((1./P_LIBOR[first_stub*sim_freq/swap.freq-1][1]-1.)/dT)
        else:
            floating.append((1./P_LIBOR[0][sim_freq/swap.freq-1]-1.)/dT)
        # then add the rest floating payments using Expectation
        for j in range(1,len(P_LIBOR[i+1]),1):
            floating.append(float((P_LIBOR[i+1][j-1]/P_LIBOR[i+1][j]-1.)/dT))
        
        if payerOrReceiver == 'payer':
            net_payment = np.asarray(floating) - c
        elif payerOrReceiver == 'receiver':
            net_payment = c - np.asarray(floating)
        else:
            print "!Error: please set the correct swap type"
    else:
        net_payment = 0  # for time maturity
    
    #first_stub = np.where(Tis>=ts[i])[0][0]
    
    # first add the first (past set) stub where T_(i-1)<t
    
    '''
    if first_stub != 0:
        floating.append((1./P_LIBOR[first_stub*sim_freq/swap.freq-1][1]-1.)/dT)
    else:
        floating.append((1./P_LIBOR[0][sim_freq/swap.freq-1]-1.)/dT)
    
    # then add the rest floating payments using Expectation
    for j in range(1,len(P_LIBOR[i+1]),1):
        floating.append(float((P_LIBOR[i+1][j-1]/P_LIBOR[i+1][j]-1.)/dT))
    '''
    '''
    if i > 118:
        print "i=",i
        print P_LIBOR[i+1]
        print first_stub
        print floating
        
    if payerOrReceiver == 'payer':
        net_payment = np.asarray(floating) - c
    elif payerOrReceiver == 'receiver':
        net_payment = c - np.asarray(floating)
    else:
        print "!Error: please set the correct swap type"
    if i == 119:
        print floating
        print P_OIS[i+1]'''
    price = swap.notional*dT*np.sum(np.asarray(P_OIS[i+1]).dot(net_payment))
    return price
    
def priceSwap_new(swap, payerOrReceiver, P_OIS, P_LIBOR, i, ts, sim_freq) :
    '''
    *** DO NOT USE THIS FUNCTION. NOT FIXED. USE priceSwap
    Args:
        swap:     a swap object
        payerOrReceiver: string type, 'payer' or 'receiver'
        P_OIS:    array, discount bond price using OIS
        P_Libor:  array, discount bond price using LIBOR
        i:        time step now, takes positive integer
        ts:       simulation time steps
        
    return:
        price:    value of swap at time ts[i]
    '''
    c = swap.coupon * swap.notional/swap.freq # fixed coupon
    dT = 1./swap.freq
    six = sim_freq/swap.freq
    first_stub = int(six-np.mod(ts[i]*sim_freq,six)-1)#####need to change
    floating = []
    # first add the first (past set) stub where T_(i-1)<t
    if first_stub != 0:
        floating.append((1./P_LIBOR[i+1+first_stub-six][six-1]-1.)/dT)
    else:
        floating.append((1./P_LIBOR[0][six-1]-1.)/dT)
    
    # then add the rest floating payments using Expectation
    for j in range(first_stub+six,len(P_LIBOR[i+1]),six):
        floating.append((P_LIBOR[i+1][j-six]/P_LIBOR[i+1][j]-1.)/dT)
        
    if payerOrReceiver == 'payer':
        net_payment = np.asarray(floating) - c
    elif payerOrReceiver == 'receiver':
        net_payment = c - np.asarray(floating)
    else:
        print "!Error: please set the correct swap type"
    
    '''firstTi = np.argmax(ts[six-1::six]>=ts[i])#####need to change this paragraph
    print "firstTi",firstTi
    print "ts[i]",ts[i]
    print P_OIS[i+1]
    startingStep = int(six-np.mod(ts[i]*sim_freq,six)-1)'''
    #print startingStep
    
    price = dT*np.sum(np.asarray(P_OIS[i+1][first_stub::six]).dot(net_payment))
    return price
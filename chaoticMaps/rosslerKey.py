#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:09:45 2022

@author: rogerchiu
"""
import numpy as np
# MaxPhase = 2.0*np.pi
MaxPhase = 2;
def rosslerKey(N,alfa = 3.8, beta = 0.05, gamma = 0.35, delta = 3.78,
                zeta = 0.2, eta = 0.1, teta = 1.9, y_1=0.3, y_2 = 0,
                y_3 = 0.05):
            
    y1 = np.empty(N*N)
    y2 = np.empty(N*N)
    y3 = np.empty(N*N)
    xy = np.empty((N,N))
    y1[0] = y_1
    y2[0] = y_2
    y3[0] = y_3
    
    for l in range((N-1)*(N-1)):
      y1[l+1] = alfa*y1[l]*(1-y1[l])-beta*(y3[l]+gamma)*(1-2*y2[l]);
      y2[l+1] = delta*y2[l]*(1-y2[l])+zeta*y3[l];
      y3[l+1] = eta*((y3[l]+gamma)*(1-2*y2[l])-1)*(1-teta*y1[l]);
    
    
    for i  in range(N-1):
        for j  in range(N-1):
            xy[i][j] = abs(y1[i*j])*MaxPhase
    return xy
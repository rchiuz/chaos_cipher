#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    chenKey(N,alpha = 1.95,beta = 1.0, x_1 = 0.3, y_1 = 0.56)
"""

import numpy as np
# MaxPhase = 2.0*np.pi
MaxPhase = 2;
def chenKey(N,alpha = 1.95,beta = 1.0, x_1 = 0.3, y_1 = 0.56 ):

#    alpha = 1.95
#    beta = 1.0
    x1 = np.empty(N*N)
    y1 = np.empty(N*N)
    xy = np.empty((N,N))
    x1[0] = x_1
    y1[0] = y_1    # con estas condiciones inicales se pude cambiar la llave
    for k in range((N-1)*(N-1)):
        x1[k+1] = 1-alpha*(x1[k]*x1[k]+y1[k]*y1[k]);
        y1[k+1] = -2*alpha*beta*x1[k]*y1[k];
    
    for i  in range(N-1):
        for j  in range(N-1):
            xy[i][j] = abs(x1[i*j])*MaxPhase
    return xy
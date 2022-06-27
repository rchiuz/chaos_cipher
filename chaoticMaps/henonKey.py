#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:11:12 2022

@author: rogerchiu
"""

import numpy as np
# MaxPhase = 2.0*np.pi
MaxPhase = 2;

def henonKey(N,aa=1.4, bb = .3, z_1 = 0.4, z_2 = 0, z_3 = 0.5):
    
    z1 = np.empty(N*N)
    z2 = np.empty(N*N)
    z3 = np.empty(N*N)
    xy = np.empty((N,N))
    z1[0]= z_1;
    z2[0]= z_2;
    z3[0]= z_3;
    
    for l in range((N-1)*(N-1)):
        z1[l+1]=aa-z2[l]*z2[l]-bb*z3[l];
        z2[l+1]=z1[l];
        z3[l+1]=z2[l];
    for i  in range(N-1):
        for j  in range(N-1):
            xy[i,j]=abs(z1[i*j])*MaxPhase
    return xy 
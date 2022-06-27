#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 09:50:23 2021

@author: rogerchiu
"""

import numpy as np
from skimage import color

def xorCipher(img, key):
    if (len(img.shape)==3):
        img = np.round(color.rgb2gray(img)*255)
    N, M = img.shape
    vectorImg = np.zeros(M*M)
    vectorKey = np.zeros(M*M)
    cipherImg = np.zeros(M*M)
    vectorImg = np.reshape(img,N*M)
    vectorKey = np.reshape(key,N*M)
    for i in range(N*M):
        cipherImg[i] = (int(vectorImg[i])^int(vectorKey[i]))
    cipherImg = np.reshape(cipherImg,[N,M])
    return cipherImg 

def xorUnCipher(cImg, key):
    N, M = cImg.shape
    vectorImg = np.zeros(M*M)
    vectorKey = np.zeros(M*M)
    unCipherImg = np.zeros(M*M)
    vectorImg = np.reshape(cImg,N*M)
    vectorKey = np.reshape(key,N*M)
    for i in range(N*M):
        unCipherImg[i] = (int(vectorImg[i])^int(vectorKey[i]))
    unCipherImg = np.reshape(unCipherImg,[N,M])
    return unCipherImg 
    
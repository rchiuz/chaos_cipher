#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from IPython import display
import config 


def resize2(img, factor):
        num = int(-np.log2(factor))
        for i in range(num):
            #img = 0.25*(img[::2,::2,...]+img[1::2,::2,...]+img[::2,1::2,...]+img[1::2,1::2,...])
            img = 0.25*(img[::2]+img[1::2])
        return img
def U_update(eta, image_est, tau):
    return SoftThresh(Psi(image_est) + eta/config.mu2, tau/config.mu2)
def SoftThresh(x, tau):
    # numpy automatically applies functions to each element of the array
    return np.sign(x)*np.maximum(0, np.abs(x) - tau)
def Psi(v):
    return np.stack((np.roll(v,1,axis=0) - v, np.roll(v, 1, axis=1) - v), axis=2)
def X_update(xi, image_est, H_fft, sensor_reading, X_divmat):
    return X_divmat * (xi + config.mu1*M(image_est, H_fft) + CT(sensor_reading))
def M(vk, H_fft):
    return np.real(fft.fftshift(fft.ifft2(fft.fft2(fft.ifftshift(vk))*H_fft)))
def C(M):
    # Image stored as matrix (row-column rather than x-y)
    top = (config.full_size[0] - config.sensor_size[0])//2
    bottom = (config.full_size[0] + config.sensor_size[0])//2
    left = (config.full_size[1] - config.sensor_size[1])//2
    right = (config.full_size[1] + config.sensor_size[1])//2
    return M[top:bottom,left:right]

def CT(b):
    v_pad = (config.full_size[0] - config.sensor_size[0])//2
    h_pad = (config.full_size[1] - config.sensor_size[1])//2
    return np.pad(b, ((v_pad, v_pad), (h_pad, h_pad)), 'constant',constant_values=(0,0))
def precompute_X_divmat(): 
    """Only call this function once! 
    Store it in a variable and only use that variable 
    during every update step"""
    return 1./(CT(np.ones(config.sensor_size)) + config.mu1)
def W_update(rho, image_est):
    return np.maximum(rho/config.mu3 + image_est, 0)
def r_calc(w, rho, u, eta, x, xi, H_fft):
    return (config.mu3*w - rho)+PsiT(config.mu2*u - eta) + MT(config.mu1*x - xi, H_fft)

def V_update(w, rho, u, eta, x, xi, H_fft, R_divmat):
    freq_space_result = R_divmat*fft.fft2( fft.ifftshift(r_calc(w, rho, u, eta, x, xi, H_fft)) )
    return np.real(fft.fftshift(fft.ifft2(freq_space_result)))
def PsiT(U):
    diff1 = np.roll(U[...,0],-1,axis=0) - U[...,0]
    diff2 = np.roll(U[...,1],-1,axis=1) - U[...,1]
    return diff1 + diff2
def MT(x, H_fft):
    x_zeroed = fft.ifftshift(x)
    return np.real(fft.fftshift(fft.ifft2(fft.fft2(x_zeroed) * np.conj(H_fft))))
def precompute_PsiTPsi():
    PsiTPsi = np.zeros(config.full_size)
    PsiTPsi[0,0] = 4 # before 4
    PsiTPsi[0,1] = PsiTPsi[1,0] = PsiTPsi[0,-1] = PsiTPsi[-1,0] = -1
    PsiTPsi = fft.fft2(PsiTPsi)
    return PsiTPsi
def precompute_R_divmat(H_fft, PsiTPsi): 
    """Only call this function once! 
    Store it in a variable and only use that variable 
    during every update step"""
    MTM_component = config.mu1*(np.abs(np.conj(H_fft)*H_fft))
    PsiTPsi_component = config.mu2*np.abs(PsiTPsi)
    id_component = config.mu3
    """This matrix is a mask in frequency space. So we will only use
    it on images that have already been transformed via an fft"""
    return 1./(MTM_component + PsiTPsi_component + id_component)
def xi_update(xi, V, H_fft, X):
    return xi + config.mu1*(M(V,H_fft) - X)

def eta_update(eta, V, U):
    return eta + config.mu2*(Psi(V) - U)

def rho_update(rho, V, W):
    return rho + config.mu3*(V - W)
def init_Matrices(H_fft):
    X = np.zeros(config.full_size)
    U = np.zeros((config.full_size[0], config.full_size[1], 2))
    V = np.zeros(config.full_size)
    W = np.zeros(config.full_size)

    xi = np.zeros_like(M(V,H_fft))
    eta = np.zeros_like(Psi(V))
    rho = np.zeros_like(W)
    return X,U,V,W,xi,eta,rho
def precompute_H_fft(psf):
    return fft.fft2(fft.ifftshift(CT(psf)))
def ADMMStep(X,U,V,W,xi,eta,rho, precomputed):
    H_fft, data, X_divmat, R_divmat = precomputed
    U = U_update(eta, V, config.tau)
    X = X_update(xi, V, H_fft, data, X_divmat)
    V = V_update(W, rho, U, eta, X, xi, H_fft, R_divmat)
    W = W_update(rho, V)
    xi = xi_update(xi, V, H_fft, X)
    eta = eta_update(eta, V, U)
    rho = rho_update(rho, V, W)
    
    return X,U,V,W,xi,eta,rho
def runADMM(psf, data, showimg = False):
    H_fft = precompute_H_fft(psf)
    X,U,V,W,xi,eta,rho = init_Matrices(H_fft)
    X_divmat = precompute_X_divmat()
    PsiTPsi = precompute_PsiTPsi()
    R_divmat = precompute_R_divmat(H_fft, PsiTPsi)
    
    for i in range(config.iters):
         X,U,V,W,xi,eta,rho = ADMMStep(X,U,V,W,xi,eta,rho, [H_fft, data, X_divmat, R_divmat])
         if i % 1 == 0:
              image = C(V)
              image[image<0] = 0
              print (i)
              if showimg:
                f = plt.figure(1)
                plt.imshow(image, cmap='gray')
                plt.title('Reconstruction after iteration {}'.format(i))
                display.display(f)
                display.clear_output(wait=True)
    return image
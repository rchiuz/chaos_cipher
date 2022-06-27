#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 11:39:43 2021
data = skimage.util.random_noise(data, mode="gaussian")
@author: rogerchiu
"""

from LightPipes import *
import matplotlib.pyplot as plt
import math
import matplotlib.image as mpimg
import numpy as np
from scipy.signal import convolve
from scipy import signal
from scipy.io import savemat
from skimage import data, color
from skimage.transform import rescale, resize, downscale_local_mean
import numpy.fft as fft
from PIL import Image
import matplotlib.pyplot as plt
from IPython import display
from skimage.color import rgb2gray
from chaoticMaps.chenKey import chenKey
from chaoticMaps.henonKey import henonKey
from chaoticMaps.rosslerKey import rosslerKey
from aDMM import *
from xOrCipher import *
import config


savedata = True #change the value of this if save results if needed


"""
Downsampling factor (used to shrink images)
"""
f = 0.25

"""
***********parameters for the optical field
"""
N=512
N2=int(N/2)
ZoomFactor=1.5
NZ=N2/ZoomFactor
a=.005*mm
wl = 532*nm
size = .02*mm
z0 = 0.1*nm;
z1 = 65*cm 
aperture = 15*mm
ap_cir = 3.*mm
seed = 7
MaxPhase = 2.0*np.pi
r = 4* mm
      
"""
"""
def mapsimage(imgDataIn, out_min,out_max):
    # rescale the final image to values between out_min,out_max in gray scale
    r_min = np.min(imgDataIn);
    r_max = np.max(imgDataIn);
    imgDataOut = np.fix((out_max-out_min)*((imgDataIn[:]-r_min)/(r_max-r_min))+(out_min));
    return imgDataOut;  
 

for alp in np.arange(1):
    alpha=2    
    """ 
    here is selected the chaotic map to be used as cipher key
    to change the map only uncoment
    """    
    key01 = chenKey(N,alpha,beta = 1.0, x_1 = 0.3, y_1 = 0.56);
    #key01 = rosslerKey(N);
    #key01 = henonKey(N,1.3);      
    
    """
          in this section the scattering medium is created using  Lightpipes library
          the phase od the medium wil be used as key
    """
      
    Fd0=Begin(35*mm,wl,N);
    Fd0 = MultPhase(Fd0,key01)
    Fd0=PipFFT(1,Fd0); # gets the fourier transform of the random field
    Fd0=CircAperture(ap_cir,0,0,Fd0); # propagates random field througt circular mask
    Fd0=PipFFT(-1,Fd0);# gets the inverse fourier transform of the  field
    I0=Intensity(0,Fd0);
    Ip0=Phase(Fd0);
    plt.imshow(I0, cmap='gray'); 
    plt.axis('off');
    plt.title('Scattering medium (difusser)')
    plt.show()
    
    
    
    Fd = RectAperture(aperture,aperture,0*mm,0*mm,0.0,Fd0)
    I0=Intensity(10,Fd);
    Ip0=Phase(Fd);
    
    
    plt.imshow(I0, cmap='gray'); 
    plt.axis('off');
    plt.title('Intensity (difusser)')
    key_2 = I0
    
    """
            here the point spread function PSF is created and 
            also is croped
    """
  
    Fps = PointSource(Fd0);
    iu = convolve(Fps.field,Ip0)
    Psf=iu[N2:(N*2)-N2,N2:(N*2)-N2]    
    plt.imshow(abs(Psf),cmap='hot'); 
    plt.axis('off');
    plt.title('Caustics at sensor plane')
    plt.show()
    Psf.shape
    
    
    """
         here the image to be ciphered is loaded and is rescaled to a range of
         -128 to 128 in gray scale
    """

    #imgname = "images/lena.png"
    #imgname = "images/cameraMan.png"
    #imgname = "images/hand.png"
    #imgname = "images/graffiti.png"
    imgname = 'images/skull2.jpg'

    
    
    img = (Image.open(imgname))
    new_image = img.resize((N, N))
    img = np.array(new_image, dtype='float32')
    img = rgb2gray(img)
    plt.imshow(img-128,cmap='gray'); plt.axis('off')
    plt.title('image to be recovery')
    plt.show()
    imgOriginal =  img-128
    
    """
        here the convolution beetwen the original image and the key
        is performed
    """

    
    iu2 = convolve(img,abs(Psf))
    imSensor=iu2[N2:(N*2)-N2,N2:(N*2)-N2]
    
    plt.imshow(imSensor,cmap='gray'); 
    plt.axis('off');
    plt.title('image at sensor plane')
    plt.show()
    imgciphered =  imSensor
    imSensor.shape
    
    imgciphered= np.uint8(mapsimage(imgciphered,0,255));
    cimg = xorCipher(imgciphered,(np.around(key01))*255)
    
    plt.imshow(cimg,cmap='gray'); 
    plt.axis('off');
    plt.title('image after XOR cipher')
    plt.show()

    
    """ here a new key will be created to be used as unciper key

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
    xy = chenKey(N,alpha,beta = 1.0, x_1 = 0.3, y_1 = 0.56); 
    
    """
          in this section the scattering medium is created using  Lightpipes library
          the phase od the medium wil be used as key
    """
    key01 = xy 
    uimg = xorUnCipher(cimg,(np.around(key01))*255)
    plt.imshow(uimg,cmap='gray'); 
    plt.axis('off');
    plt.title('image after XOR uncipher')
    plt.show()
    imSensor2 = uimg
    imSensor2 = imSensor
      
    Fd0=Begin(35*mm,wl,N);
    Fd0 = MultPhase(Fd0,xy)
    Fd0=PipFFT(1,Fd0); # gets the fourier transform of the random field
    Fd0=CircAperture(ap_cir,0,0,Fd0); # propagates random field througt circular mask
    Fd0=PipFFT(-1,Fd0);# gets the inverse fourier transform of the  field
    I0=Intensity(0,Fd0);
    Fd = RectAperture(aperture,aperture,0*mm,0*mm,0.0,Fd0)
    I0=Intensity(10,Fd);
    Ip0=Phase(Fd);
    
    
    Fps = PointSource(Fd0);
    iu = convolve(Fps.field,Ip0)
    Psf=iu[N2:(N*2)-N2,N2:(N*2)-N2]
    Psf=np.array(Psf)
    dataSensor=np.array(imSensor2)
    psf = np.array(Psf.real, dtype='float32')
    data = np.array(dataSensor.real, dtype='float32')
    
    bg = np.mean(psf[5:6,5:6]) 
    bg = .001
    psf -= bg
    data -= bg
    psf /= np.linalg.norm(psf.ravel())
    data /= np.linalg.norm(data.ravel())
    config.sensor_size = np.array(psf.shape)
    config.full_size = 2*config.sensor_size
    
    
    """
        here the ADMM algorithm is used to uncipher the image
    """
    config.iters = 200
    final_im = runADMM(psf, data, showimg = True) #change showimg = False to avoid showing image of each itaration    
    nfinal_im = np.uint8(mapsimage(final_im,0,255));
    
    plt.imshow(nfinal_im, cmap='gray')
    plt.title('Final reconstructed image after {} iterations'.format(config.iters))
    display.display()
    plt.axis('off');
    
    
    imgOriginal = np.uint8(mapsimage(imgOriginal,0,255));
    imgciphered= np.uint8(mapsimage(imgciphered,0,255));
    
    
    
    """
        here the output varibles are recorded in .mat format 
    """
    if savedata:
       imgoutName =  imgname.split('.')[0]
       mdic = {"key":psf,"imgOriginal":imgOriginal , "imgCiphered":imgciphered,"imgXORciphered":cimg,"imgRetrieved":nfinal_im }
       filename = "outputdata/{0}_exp_{1}.mat".format(imgoutName.split('/')[1],alpha)
       savemat(filename, mdic)



    
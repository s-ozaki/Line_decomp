from astropy.io import fits
import numpy as np
import scipy.constants as scon

def GetLambdaArr(hdr):
    # Creating array of lambdas
    # hdl: FITS hdl
    
    #npix = hdl[0].data.shape[0]
    if 'DISPAXIS' in hdr:
        dispax = hdr['DISPAXIS']
    else:
        dispax = 3
        
    npix = hdr['NAXIS'+str(dispax)]
    pixarr = np.array(range(npix), dtype=np.int) + 1
    lam = Pix2Wav(pixarr, hdr)
    return lam
    
def Pix2Wav(pix, hdr):
    if 'DISPAXIS' in hdr:
        dispax = hdr['DISPAXIS']
    else:
        dispax = 3
    crval = hdr['CRVAL'+str(dispax)]
    crpix = hdr['CRPIX'+str(dispax)]
    cdelt_key = 'CD'+str(dispax)+'_'+str(dispax)
    if cdelt_key in hdr:
        cdelt = hdr[cdelt_key]
    else:
        cdelt = hdr['CDELT'+str(dispax)]
    
    w = crval + (pix - crpix + 1) * cdelt
    return w


def Wav2Pix(w, hdr):
    # w: wavelength (A)
    
    if 'DISPAXIS' in hdr:
        dispax = hdr['DISPAXIS']
    else:
        dispax = 3
    crval = hdr['CRVAL'+str(dispax)]
    crpix = hdr['CRPIX'+str(dispax)]
    cdelt_key = 'CD'+str(dispax)+'_'+str(dispax)
    if cdelt_key in hdr:
        cdelt = hdr[cdelt_key]
    else:
        cdelt = hdr['CDELT'+str(dispax)]
    
    pix = (w - crval) / cdelt + crpix - 1.0
    return pix


def Wav2Vel(w, rw, z):
    # w: wavelength (A)
    # rw: rest wavelength (A)
    # z: redshift
    
    v = (w / (1.0 + z) / rw - 1.0)*scon.c/1000.0
    return v  # (km/s)


def Vel2Wav(v, rw, z):
    # v: velocity (km/s)
    # rw: rest wavelength (A)
    # z: redshift
    
    w = (1.0 + v*1000.0 / scon.c) * (1.0 + z) * rw
    return w  # (A)
    


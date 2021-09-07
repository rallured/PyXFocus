import numpy as np
import matplotlib.pyplot as plt
import reconstruct
import pdb  

def padArrays(imglist):
    """
    Pad arrays with a border of 100s, which the reconstructor
    function interprets as NaNs
    """
    outimg = []
    for i in imglist:
        sh = np.shape(i)
        temp = np.zeros((sh[0]+2,sh[1]+2),order='F')+100.
        temp[1:-1,1:-1] = i
        outimg.append(temp)
        
    return outimg

def southwell(gx,gy,criteria,h,maxiter=10000):
    """
    Take 2D x and y gradient arrays and perform Southwell reconstruction.
    If RMS difference in slope from one iteration to the next falls below
    criteria, then algorithm stops
    gx and gy should be in units of slope (rise/run)
    h is in units of step size (e.g. lenslet width)
    """
    #Handle missing data
    ind = np.logical_or(np.isnan(gx),np.isnan(gy))
    gx[ind] = 100.
    gy[ind] = 100.
    #Create phase object
    phase = np.zeros(np.shape(gx),order='F')
    phase[ind] = 100.
    #Pad arrays
    phase,gx,gy = padArrays([phase,gx,gy])
    #Send through reconstructor
    reconstruct.reconstruct(gx,gy,1e-10,1.,phase,maxiter)
    #Strip off border
    phase = phase[1:-1,1:-1]
    #Handle missing data
    phase[ind] = np.nan
    #Return inverse due to reconstructor sign convention
    return -phase
    

def example():
    """
    Example of how the reconstructor works
    """
    #Create test image
    xg,yg = np.meshgrid(np.linspace(-1,1,100),np.linspace(-1,1,100))
    img = np.polynomial.legendre.legval2d(xg,yg,[[0,1,0],[0,.5,0],[1,0,0]])
    #Take gradients
    gx,gy = np.gradient(img)
    #Create a boundary region
    rad = np.sqrt(xg**2+yg**2)
    img[rad>1] = np.nan
    gx[rad>1] = np.nan
    gy[rad>1] = np.nan
    #Reconstruct wavefront
    recon = southwell(gx,gy,1e-12,1.)

    #Plot results
    fig = plt.figure()
    fig.add_subplot(1,3,1)
    plt.imshow(img)
    plt.title('Original')
    plt.colorbar()
    fig.add_subplot(1,3,2)
    plt.imshow(recon)
    plt.title('Reconstructed')
    plt.colorbar()
    fig.add_subplot(1,3,3)
    plt.imshow(img-recon)
    plt.title('Residual')
    plt.colorbar()
    return recon

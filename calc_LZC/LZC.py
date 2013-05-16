import numpy as np
import pickle,os


## Load diagnostic calibration data
#try:
     #import pkgutil
     #data = pkgutil.get_data(__name__, 'LZCpar.p')
#except ImportError:
     #import pkg_resources
     #data = pkg_resources.resource_string(__name__, 'LZCpar.p')
#LZCp=pickle.load(data)

#LZCp=pickle.load(open(os.path.dirname(__file__)+'/LZCpar.p' ,'r'))
LZCp=np.load(os.path.dirname(__file__)+'/LZC2_all.npy')

def cLZC(u=None,g=None,r=None,i=None,z=None,du=None,dg=None,dr=None,di=None,dz=None,M='g',col='r',diag='T04',pmeth='model',N=1000):
    """
    Estimate galaxy metallicity based on two bands of galaxy photometry using 
    the LZC relation
    
    If (dg,dr,di) photometric uncertainties are supplied, a Monte Carlo calculation 
    will be performed and uncertainties will be returned
    
    Keyword arguments:
    u,g,r,i,z -- K-corrected galaxy absolute magnitudes in Sloan filters (floats or 1D arrays)
    du,dg,dr,di,dz -- Corresponding photometric uncertainties  (floats or 1D arrays)
    M -- ugriz band to use as luminosity (default 'g')
    col -- ugriz band to use for color (default 'r', meaning g-r - must be redder than M filter)
    diag -- Metallicity diagnostic to return result in (default 'T04', options below)
    pmeth -- Photometric system (default 'model'. options below)
    N -- Number of Monte Carlo samples to use for uncertainty calculation (default 1000)
    
    Available metallicity diagnostics:
    * M91
    * P01 (poorly calibrated)
    * KD02
    * KK04
    * PP04_N2
    * PP04_O3N2
    * T04
    * PT05 (poorly calibrated)
    * PVT10
    
    Available photometric systems:
    * petro (Petrosian)
    * model
    * cmodel
    * fiber (3 arcsecond - u and z bands not avaialable)
    
    Returns: (Z,dZ,disp)
    * Z: Point estimate (median) for metallicity in the specified diagnostic (T04 by default).
    * dZ: Metallicity uncertainty according to propogation of photometric uncertainty 
        (np.nan if dg,dr,di not specified).    Think of this as a statistical uncertainty.
    * disp: Intrinsic scatter (standard deviation of Metallicity residuals) for this 
        calibration of the LZC.  Think of this as a systematic uncertainty.
    * warnings: A binary array indicating if (0) the galaxy is in the calibrated mu range
       or (1) the galaxy is outside the calibrated range and the LZC was extrapolated
    
    For more information: https://www.cfa.harvard.edu/~nsanders/papers/LZC/summary.htm
    """
    ## Load calibration data
    #balpha,disp0,dispb,p,disp,mu_min,mu_max=LZCp[diag][M][col]
    sel=np.where((LZCp['diag']==diag) & (LZCp['pmeth']==pmeth) & (LZCp['f1']==M) & (LZCp['f2']==col))
    balpha=LZCp[sel]['balpha']
    mu_min=LZCp[sel]['mu_min']
    mu_max=LZCp[sel]['mu_max']
    dispb=LZCp[sel]['dispb']
    p=[LZCp[sel]['p0'],LZCp[sel]['p1'],LZCp[sel]['p2'],LZCp[sel]['p3']]
    ## Load input data
    phot={'u':u,'g':g,'r':r,'i':i,'z':z,'du':du,'dg':dg,'dr':dr,'di':di,'dz':dz}
    # If inputs aren't lists, convert them
    for key in phot:
        try:
            phot[key][0]
        except (TypeError,IndexError) as e:
            phot[key]=np.array([phot[key]])
    ## Establish Monte Carlo samples
    MCphot_1=np.ones([len(phot[M]),N])
    MCphot_2=np.ones([len(phot[M]),N])
    MCphot_3=np.ones([len(phot[M]),N])
    for i in range(len(phot[M])):
        if phot['d'+M]==[None]: 
            MCphot_1[i]*=phot[M][i]
            MCphot_2[i]*=phot[M][i]
            MCphot_3[i]*=phot[col][i]
        else: 
            MCphot_1[i]=np.random.normal(phot[M][i],phot['d'+M][i],N)
            MCphot_2[i]=np.random.normal(phot[M][i],phot['d'+M][i],N)
            MCphot_3[i]=np.random.normal(phot[col][i],phot['d'+col][i],N)
    ## Apply Z calibration
    cZ=lambda x: np.polyval(p,x)
    mu=MCphot_1-balpha*(MCphot_2-MCphot_3)
    MC_Z=cZ(mu)
    Z=np.median(MC_Z,axis=1)
    dZ=np.std(MC_Z,axis=1)
    ## Replace entries without uncertainty with np.nan
    sel=np.where(dZ<1e-7)
    dZ[sel]=[np.nan]*len(sel[0])
    ##Issue warnings if outside of calibrated range
    warn=(((np.median(mu,axis=1)>mu_min) & (np.median(mu,axis=1)<mu_max))==False)*1
    #if len(warn)>0: warnings.warn('The following objects are outside the calibrated LZC range:\n'+'\n'.join(warn.astype('str')))
    return Z,dZ,dispb,warn
    

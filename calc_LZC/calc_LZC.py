import numpy as np
import pickle,warnings

## Load diagnostic calibration data
LZCpar=pickle.load(open(LZCpar.p,'r'))

def cLZC(g=None,r=None,i=None,dg=None,dr=None,di=None,M='g',col='gr',diag='T04',N=1000):
    """
    Estimate galaxy metallicity based on two bands of galaxy photometry using 
    the LZC relation
    
    If (dg,dr,di) photometric uncertainties are supplied, a Monte Carlo calculation 
    will be performed and uncertainties will be returned
    
    Keyword arguments:
    g,r,i -- K-corrected galaxy absolute magnitudes in Sloan filters (1D arrays)
    dg,dr,di -- Corresponding photometric uncertainties  (1D arrays)
    M -- gri band to use as luminosity (default 'g')
    col -- gri color to use (default 'gr', meaning g-r)
    diag -- Metallicity diagnostic to return result in (default 'T04', options below)
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
    
    Returns: (Z,dZ,disp)
    * Z: Point estimate (median) for metallicity in the specified diagnostic (T04 by default).
    * dZ: Metallicity uncertainty according to propogation of photometric uncertainty 
        (nan if dg,dr,di not specified).    Think of this as a statistical uncertainty.
    * disp: Intrinsic scatter (standard deviation of Metallicity residuals) for this 
        calibration of the LZC.  Think of this as a systematic uncertainty.
    
    For more information: https://www.cfa.harvard.edu/~nsanders/papers/LZC/summary.htm
    """
    ## Load calibration data
    balpha,disp0,dispb,p,disp,mu_min,mu_max=LZCp[diag][M][col]
    ## Load input data
    phot={'g':g,'r':r,'i':i,'dg':dg,'dr':dr,'di':di}
    ## Establish Monte Carlo samples
    MCphot_1=ones([len(phot[M]),N])
    MCphot_2=ones([len(phot[M]),N])
    MCphot_3=ones([len(phot[M]),N])
    for i in range(len(phot[M])):
        if phot['d'+M][i]==None: 
            MCphot_1[i]*=phot['d'+M][i]
            MCphot_2[i]*=phot['d'+col[0]][i]
            MCphot_3[i]*=phot['d'+col[1]][i]
        else: 
            MCphot_1[i]=np.random.normal(phot[M][i],phot['d'+M][i],N)
            MCphot_2[i]=np.random.normal(phot[col[0]][i],phot['d'+col[0]][i],N)
            MCphot_3[i]=np.random.normal(phot[col[1]][i],phot['d'+col[1]][i],N)
    ## Apply Z calibration
    cZ=lambda x: polyval(p,x)
    mu=MCphot_1-balpha*(MCphot_2-MCphot_3)
    MC_Z=cZ(mu)
    Z=median(MC_Z,axis=1)
    dZ=std(MC_Z,axis=1)
    ## Replace entries without uncertainty with nan
    sel=where(dZ==0)
    dZ[sel]=[nan]*len(sel[0])
    ##Issue warnings if outside of calibrated range
    warn=where(((median(mu,axis=1)>mu_min) & (median(mu,axis=1)<mu_max))==False)[0]
    if len(warn)>0: warnings.warn('The following objects are outside the calibrated LZC range:\n'+'\n'.join(warn.astype('str')))
    return Z,dZ,dispb
    

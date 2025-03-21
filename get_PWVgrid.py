## This script is a function / definition of calculating PWV grids
## Created by C. Bayu Risanto, S.J. (27 November 2023)
import numpy as np
import random
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata
from netCDF4 import Dataset,num2date
import time
import glob

##### define the gridded PWV (whole domain) Function ###########################################
def get_griddedPWV(ncfile):
    dx = xr.open_dataset(ncfile)
    QVAPOR = dx['QVAPOR'].squeeze()
    T = dx['T']
    PH = dx['PH'].squeeze()
    PHB = dx['PHB'].squeeze()
    MU = dx['MU']
    MUB = dx['MUB']
    DNW = dx['DNW'].squeeze()
    Q2 = dx['Q2']
    PSFC = dx['PSFC']
    XLON = dx.XLONG
    XLAT = dx.XLAT
    nk = len(DNW)
    nt,ny,nx = MU.shape

## calculate total pressure on mass point (half(mass) levels, T-point)
    gas_constant   = 287.0
    gas_constant_v = 461.6
    ps0            = 100000.0
    ts0            = 300.0
    
    rd_over_rv = gas_constant / gas_constant_v
    cpovcv = 1.4

#Adapted the code from WRF module_big_step_utilities_em.F ----
#         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)
# Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

    qvf1 = 1.0 + QVAPOR / rd_over_rv
    rhoo = [ -(MUB+MU)/ (( (PH[k+1,:,:] + PHB[k+1,:,:]) - (PH[k,:,:] + PHB[k,:,:]) ) / DNW[k]) for k in range(nk) ]
    rho = np.squeeze(np.asarray(rhoo))

# .. total pressure
    totalPress = ps0 * ( (gas_constant*(ts0+T)*qvf1) / (ps0/rho) )**cpovcv
    totalPs = np.squeeze(totalPress)

# get column of pressure and specific humidity
    Pall = np.concatenate((PSFC,totalPs),axis=0)
    QVall = np.concatenate((Q2 / (1.0 + Q2),QVAPOR / (1.0 + QVAPOR)),axis=0)

# get PWV similar in DART
    PWV1 = [[[ 0.5*(QVall[k,j,i] + QVall[k+1,j,i]) * (Pall[k,j,i]-Pall[k+1,j,i]) 
            for i in range(nx)] for j in range(ny)] for k in range(nk)]
    PWV = np.sum(np.asarray(PWV1),axis=0)/9.81
    return PWV,XLON,XLAT



##### modified with absolute values #######################################################################

def get_griddedPWV_abs(ncfile):
    dx = xr.open_dataset(ncfile)
    QVAPOR = dx['QVAPOR'].squeeze()
    T = dx['T']
    PH = dx['PH'].squeeze()
    PHB = dx['PHB'].squeeze()
    MU = dx['MU']
    MUB = dx['MUB']
    DNW = dx['DNW'].squeeze()
    Q2 = dx['Q2']
    PSFC = dx['PSFC']
    XLON = dx.XLONG
    XLAT = dx.XLAT
    nk = len(DNW)
    nt,ny,nx = MU.shape

## calculate total pressure on mass point (half(mass) levels, T-point)
    gas_constant   = 287.0
    gas_constant_v = 461.6
    ps0            = 100000.0
    ts0            = 300.0
    
    rd_over_rv = gas_constant / gas_constant_v
    cpovcv = 1.4

#Adapted the code from WRF module_big_step_utilities_em.F ----
#         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)
# Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)
# /st8/arellano/DARTM/DART/observations/forward_operators/obs_def_tpw_mod.f90
# Modified by C. Bayu Risanto, S.J & A.F.Arellano

    qvf1 = 1.0 + QVAPOR / rd_over_rv
    rhoo = [ -(MUB+MU)/ (( (PH[k+1,:,:] + PHB[k+1,:,:]) - (PH[k,:,:] + PHB[k,:,:]) ) / DNW[k]) for k in range(nk) ]
    rho = np.squeeze(np.asarray(rhoo))

# .. total pressure
    totalPress = ps0 * ( (gas_constant*(ts0+T)*qvf1) / (ps0/rho) )**cpovcv
    totalPs = np.squeeze(totalPress)

# get column of pressure and specific humidity
    Pall = np.concatenate((PSFC,totalPs),axis=0)
    QVall = np.concatenate((Q2 / (1.0 + Q2),QVAPOR / (1.0 + QVAPOR)),axis=0)

# get PWV similar in DART
    PWV1 = [[[ 0.5*np.absolute(QVall[k,j,i] + QVall[k+1,j,i]) * np.absolute(Pall[k,j,i]-Pall[k+1,j,i]) 
            for i in range(nx)] for j in range(ny)] for k in range(nk)]
    PWV = np.sum(np.asarray(PWV1),axis=0)/9.81
    return PWV,XLON,XLAT

## This following functions are to calculate statical skills of RAINNC
## i.e., CSI, POD, FAR
## Created by C. Bayu Risanto, S.J. (25 June 2024)

import numpy as np

## 3 by 3 neighborhood evaluation ###########################################
def F_rainmetrics_NV9(m_RAIN_X,m_OBSV_X,threshold):
    shape = m_RAIN_X.shape
    ntime = shape[0]; nlat = shape[1]; nlon = shape[2]
    #print(ntime);print(nlat);print(nlon)

    stat = np.zeros(shape=(ntime,4,nlat,nlon))  # create multi-dim array

    for time in range(ntime):
        #print(time)
        for j in range(1,nlat-1):
            for i in range(1,nlon-1):
                count1 = 0
                count2 = 0

                for jc in range(-1,2):
                    for ic in range(-1,2):
                        if m_RAIN_X[time,j+jc,i+ic] >= threshold:
                            count1 += 1
                        if m_OBSV_X[time,j+jc,i+ic] >= threshold:
                            count2 += 1

                if count1 > 0:
                    mod = True
                else:
                    mod = False
                if count2 > 0:
                    obs = True
                else:
                    obs = False

                # calculate hits, misses, false alarms, and corr negatives
                if obs == True and mod == True:
                    #hits
                    stat[time,0,j,i] = stat[time,0,j,i] + 1
                elif obs == True and mod == False:
                    #miss
                    stat[time,1,j,i] = stat[time,1,j,i] + 1
                elif obs == False and mod == True:
                    #FA
                    stat[time,2,j,i] = stat[time,2,j,i] + 1
                elif obs == False and mod == False:
                    #CN
                    stat[time,3,j,i] = stat[time,3,j,i] + 1

    # summary
    hits = np.nansum(stat[:,0,:,:],axis=0)
    miss = np.nansum(stat[:,1,:,:],axis=0)
    fals = np.nansum(stat[:,2,:,:],axis=0)
    corn = np.nansum(stat[:,3,:,:],axis=0)

    # metrics
    CSI = hits/(hits+fals+miss)
    POD = hits/(hits+miss)
    FAR = fals/(hits+fals)
    POFD = fals/(fals+corn)
    return CSI, POD, FAR, stat

## 5 by 5 neighborhood evaluation ############################################
def F_rainmetrics_NV25(m_RAIN_X,m_OBSV_X,threshold):
    shape = m_RAIN_X.shape
    ntime = shape[0]; nlat = shape[1]; nlon = shape[2]
    #print(ntime);print(nlat);print(nlon)

    stat = np.zeros(shape=(ntime,4,nlat,nlon))  # create multi-dim array

    for time in range(ntime):
        #print(time)
        for j in range(2,nlat-2):
            for i in range(2,nlon-2):
                count1 = 0
                count2 = 0

                for jc in range(-2,3):
                    for ic in range(-2,3):
                        if m_RAIN_X[time,j+jc,i+ic] >= threshold:
                            count1 += 1
                        if m_OBSV_X[time,j+jc,i+ic] >= threshold:
                            count2 += 1

                if count1 > 0:
                    mod = True
                else:
                    mod = False
                if count2 > 0:
                    obs = True
                else:
                    obs = False

                # calculate hits, misses, false alarms, and corr negatives
                if obs == True and mod == True:
                    #hits
                    stat[time,0,j,i] = stat[time,0,j,i] + 1
                elif obs == True and mod == False:
                    #miss
                    stat[time,1,j,i] = stat[time,1,j,i] + 1
                elif obs == False and mod == True:
                    #FA
                    stat[time,2,j,i] = stat[time,2,j,i] + 1
                elif obs == False and mod == False:
                    #CN
                    stat[time,3,j,i] = stat[time,3,j,i] + 1

    # summary
    hits = np.nansum(stat[:,0,:,:],axis=0)
    miss = np.nansum(stat[:,1,:,:],axis=0)
    fals = np.nansum(stat[:,2,:,:],axis=0)
    corn = np.nansum(stat[:,3,:,:],axis=0)

    # metrics
    CSI = hits/(hits+fals+miss)
    POD = hits/(hits+miss)
    FAR = fals/(hits+fals)
    POFD = fals/(fals+corn)
    return CSI, POD, FAR, stat

## 7 by 7 neighborhood evaluation #########################################
def F_rainmetrics_NV49(m_RAIN_X,m_OBSV_X,threshold):
    shape = m_RAIN_X.shape
    ntime = shape[0]; nlat = shape[1]; nlon = shape[2]
    #print(ntime);print(nlat);print(nlon)
    
    stat = np.zeros(shape=(ntime,4,nlat,nlon))  # create multi-dim array
    
    for time in range(ntime):
        #print(time)
        for j in range(3,nlat-3):
            for i in range(3,nlon-3):
                count1 = 0
                count2 = 0
            
                for jc in range(-3,4):
                    for ic in range(-3,4):
                        if m_RAIN_X[time,j+jc,i+ic] >= threshold:
                            count1 += 1
                        if m_OBSV_X[time,j+jc,i+ic] >= threshold:
                            count2 += 1
                
                if count1 > 0:
                    mod = True
                else:
                    mod = False
                if count2 > 0:
                    obs = True
                else:
                    obs = False
                
                # calculate hits, misses, false alarms, and corr negatives
                if obs == True and mod == True:
                    #hits
                    stat[time,0,j,i] = stat[time,0,j,i] + 1 
                elif obs == True and mod == False:
                    #miss
                    stat[time,1,j,i] = stat[time,1,j,i] + 1
                elif obs == False and mod == True:
                    #FA
                    stat[time,2,j,i] = stat[time,2,j,i] + 1
                elif obs == False and mod == False:
                    #CN
                    stat[time,3,j,i] = stat[time,3,j,i] + 1

    # summary
    hits = np.nansum(stat[:,0,:,:],axis=0)
    miss = np.nansum(stat[:,1,:,:],axis=0)
    fals = np.nansum(stat[:,2,:,:],axis=0)
    corn = np.nansum(stat[:,3,:,:],axis=0)
    
    # metrics
    CSI = hits/(hits+fals+miss)
    POD = hits/(hits+miss)
    FAR = fals/(hits+fals)
    POFD = fals/(fals+corn)
    return CSI, POD, FAR, stat

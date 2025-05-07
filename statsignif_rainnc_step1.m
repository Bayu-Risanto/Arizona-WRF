close all; clear all; clc;
% This script is to randomize the DA and noDA rainnc and calculate the CSI,
% POD, and FAR difference. Then save them for field significant analyses.
% Created  by C. Bayu Risanto, S.J. (26 June 2024)
%% read data
link = '/br1/castrogroup/bayu/UA-WRF/ARIZONA_PRO/WRF_EXTR/regRAINNC/summary/';
f18d = load([link '18hour_accum_assim.RAINNC+OBS_AZ+.mat']);
f18n = load([link '18hour_accum_noassim.RAINNC+OBS_AZ+.mat']);
f03d = load([link '3hour_accum_assim.RAINNC+OBS_AZ+.mat']);
f03n = load([link '3hour_accum_noassim.RAINNC+OBS_AZ+.mat']);

r18d = double(f18d.RAINNC18h);
r18n = double(f18n.RAINNC18h);
o18h = double(f18d.Obs18h);

r3d = double(f03d.RAINNC0912);  %% change here 1821,2100,0003,0306,0609,0912
r3n = double(f03n.RAINNC0912);  %% change here too
o3h_ALL = double(f03d.Obs3hrly); 
o3h = squeeze(o3h_ALL(6,:,:,:)); %% change here too 1,2,3,4,5,6

% merge and construct   changer here!!!!!!!!!!!!!!!!!!!!!!!!
WRF = cat(1,r3d,r3n);
OBS = repmat(o3h,2,1,1);

% threshold    change here !!!!!!!!!!!!!!!!!!!!!!!
threshold_st = 1.9;

% file out     change here !!!!!!!!!!!!!!!!!!!!
fout = 'diff_0912_total_1000CSIPODFAR_NV25_rainfall_thr';

%% MONTE CARLO here
x = 1000;
perts = NaN(x,24+24);
rng('default')
for i = 1:x
    perts(i,:) = randperm(48); 
end

threshold = threshold_st;
thresh = sprintf('%3.1f',threshold);

nlon = length(r18d(1,1,:));
nlat = length(r18d(1,:,1));

diff_CSI = NaN(nlon,nlat,x);
diff_POD = diff_CSI;
diff_FAR = diff_CSI;

for i = 1:x
    
   
    WRF_wa = WRF(perts(i,1:24),:,:);
    WRF_wn = WRF(perts(i,25:end),:,:);
    OBS_wa = OBS(perts(i,1:24),:,:);
    OBS_wn = OBS(perts(i,25:end),:,:);
    
    % analysis
    
%     [CSI_wa,POD_wa,FAR_wa] = F_calc_rainmetrics_24h(WRF_wa,GPM_wa,threshold); %This is for per one grid evaluation 
%     [CSI_wn,POD_wn,FAR_wn] = F_calc_rainmetrics_24h(WRF_wn,GPM_wn,threshold); %This is for per one grid evaluation 
    
    [CSI_wa,POD_wa,FAR_wa] = F_rainmetrics_NV25(WRF_wa,OBS_wa,threshold);
    [CSI_wn,POD_wn,FAR_wn] = F_rainmetrics_NV25(WRF_wn,OBS_wn,threshold);
    
    diff_CSI(:,:,i) = CSI_wa - CSI_wn;
    diff_POD(:,:,i) = POD_wa - POD_wn;
    diff_FAR(:,:,i) = FAR_wa - FAR_wn;
    
   count = i
end

%% LET's save this 1000 randomly picked strong and weak case difference
dir_out = '/br1/castrogroup/bayu/UA-WRF/ARIZONA_PRO/WRF_EXTR/CSIPODFAR1000_AZ+/'; %Change here!!!!
save ([dir_out fout thresh '.mat'], 'diff_CSI', 'diff_POD','diff_FAR','-v7.3'); 


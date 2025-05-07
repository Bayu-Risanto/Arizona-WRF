close all;clear all;clc;
% This script is to analyze the CSI POD FAR using statistical field significance
% by taking the layers above 95% and below 50 by Livezey and Chen (1983) and Monte Carlo technique.
% input is from statsignif_rainnc_step1.m
% output is the percentage of field significant.
%% read data
tp = 7; %Change here!!!! 1 to 7
group = {'1821';'2100';'0003';'0306';'0609';'0912';'18'}; 
nv = '25';
threshold_st = [2.2; 3.6; 3.4; 2.9; 2.6; 1.9; 9.4];
threshold = threshold_st(tp);
thresh = sprintf('%3.1f',threshold);

%% load the 1000 diff_CSIPODFAR
dir_in = '/br1/castrogroup/bayu/UA-WRF/ARIZONA_PRO/WRF_EXTR/CSIPODFAR1000_AZ+/';
load ([dir_in 'diff_',group{tp},'_total_1000CSIPODFAR_NV' nv '_rainfall_thr' thresh '.mat']);

%% SORT all of the pixel; each pixel has 1000 diff variance
sortCSI = sort(diff_CSI,3);
sortPOD = sort(diff_POD,3);
sortFAR = sort(diff_FAR,3);

%% get the 50th and 950th values from each pixel for NV25
tail_CSI_1 = sortCSI(:,:,50); %50
tail_CSI_2 = sortCSI(:,:,950); %950
tail_POD_1 = sortPOD(:,:,50);
tail_POD_2 = sortPOD(:,:,950);
tail_FAR_1 = sortFAR(:,:,1);
tail_FAR_2 = sortFAR(:,:,900);

%% then compare to the original difference; take the sum
dir_ORI = '/br1/castrogroup/bayu/UA-WRF/ARIZONA_PRO/WRF_EXTR/CSIPODFAR_AZ+/';
difORI = load ([dir_ORI 'diff_',group{tp},'_CSIPODFAR_NV' nv '.rainfall_thr' thresh '.mat']);
csi_val = permute(difORI.DIFF_CSI,[2,1]);
pod_val = permute(difORI.DIFF_POD,[2,1]);
far_val = permute(difORI.DIFF_FAR,[2,1]);

nlat = length(diff_CSI(:,1,1));
nlon = length(diff_CSI(1,:,1));

for i = 1:nlat
    for j = 1:nlon
        if csi_val(i,j) >= tail_CSI_2(i,j) || csi_val(i,j) <= tail_CSI_1(i,j)
            test_csi(i,j) = 1;
        else
            test_csi(i,j) = 0;
        end
        
        if pod_val(i,j) >= tail_POD_2(i,j) || pod_val(i,j) <= tail_POD_1(i,j)
            test_pod(i,j) = 1;
        else
            test_pod(i,j) = 0;
        end
        
        if far_val(i,j) > tail_FAR_2(i,j) || far_val(i,j) < tail_FAR_1(i,j)
            test_far(i,j) = 1;
        else
            test_far(i,j) = 0;
        end
        
    end
    
end
origmap_csi = nansum(nansum(test_csi));
origmap_pod = nansum(nansum(test_pod));
origmap_far = nansum(nansum(test_far));

%% then compare to each layer in the 1000 differnce; take the sum 
for k = 1:length(sortCSI(1,1,:))
    
    for i = 1:nlat
        for j = 1:nlon

            if diff_CSI(i,j,k) >= tail_CSI_2(i,j) || diff_CSI(i,j,k) <= tail_CSI_1(i,j)
                test_csi_1000(i,j) = 1;
            else
                test_csi_1000(i,j) = 0;
            end
            
            if diff_POD(i,j,k) >= tail_POD_2(i,j) || diff_POD(i,j,k) <= tail_POD_1(i,j)
                test_pod_1000(i,j) = 1;
            else
                test_pod_1000(i,j) = 0;
            end
            
            if diff_FAR(i,j,k) > tail_FAR_2(i,j) || diff_FAR(i,j,k) < tail_FAR_1(i,j)
                test_far_1000(i,j) = 1;
            else
                test_far_1000(i,j) = 0;
            end
            
        end
        
    end
    
    total_csi(k) = nansum(nansum(test_csi_1000));
    total_pod(k) = nansum(nansum(test_pod_1000));
    total_far(k) = nansum(nansum(test_far_1000));
    disp(k)
end

%% sort the 1000 values
sorted_1000csi = (sort(total_csi))';
sorted_1000pod = (sort(total_pod))';
sorted_1000far = (sort(total_far))';

%% get the percentage

sorted = horzcat(sorted_1000csi, sorted_1000pod, sorted_1000far); 
orig   = [origmap_csi; origmap_pod; origmap_far];

for i = 1:3
    
    IDX = find( (abs(sorted(:,i) - orig(i)) == min(abs(sorted(:,i) - orig(i))) ) );
    % get the higher IDX here
    idx = IDX(end,:);

    expan = 1;  % change here if necessary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    x = idx-expan:idx+expan;
    
    if idx == 1000
        signif(i) = 99.9;
    elseif idx == 1
           signif(i) = 0;
    else
        signif(i) = idx/10;
    end
        
%     elseif sorted(idx,i) == orig(i)
%            signif(i) = round( idx/10 , 1);
%     elseif sorted(idx,i) > orig(i)
%             num = orig(i) - sorted(idx-expan,i) ;
%             denum = sorted(idx,i) - sorted(idx-expan,i);
%             signif(i) = round( ((idx-expan) + (num/denum))/10 , 1);
%     elseif sorted(idx,i) < orig(i)
%             num = orig(i) - sorted(idx,i);
%             denum = sorted(idx+expan,i) - sorted(idx,i);
%             signif(i) = round ( ((idx) + (num/denum))/10 , 1); 
%     end
    
end

disp(signif)

%% save the significance field percentage
dir_out = '/br1/castrogroup/bayu/UA-WRF/ARIZONA_PRO/WRF_EXTR/FIELDSIGNIF_AZ+/';
save ([dir_out 'signif_field_' group{tp} '_rainfall_' thresh '.mat'], 'signif','test_csi','test_pod','test_far');


function [CSI,POD,FAR] = F_rainmetrics_NV25(WRF,observed,threshold)
[ntime,nlat,nlon] = size(observed);

% This fuction is for 25 grids neighborhood verification
stats1 = false(nlon, nlat, 4, ntime);

for time = 1:ntime
    
    for i = 3:nlon-2
        %disp(i);
        for j = 3:nlat-2
            if observed(time,j,i) >= threshold
                obs = true;
            else
                obs = false;
            end
            
            %start looping through the 25 surrounding cells
            count1 = 0;
            count2 = 0;
            
            for ic = -2:2
                for jc = -2:2
                    if WRF(time,j+jc,i+ic) >= threshold
                        count1 = count1 + 1;
                    end
                     if observed(time,j+jc,i+ic) >= threshold
                        count2 = count2 + 1;
                    end
                end
            end
            
            if count1 > 0
                mod = true;
            else
                mod = false;
            end
        
            if count2 > 0
                obs2 = true;
            else
                obs2 = false;
            end
            
            % now calculate hits, misses, false alarm, and correct negative
            if obs2 == true && mod == true
                %hit
                stats1(i,j,1,time) = true;
            elseif obs2 == true && mod == false
                % miss
                stats1(i,j,2,time) = true;
            elseif obs2 == false && mod == true
                % false alarm
                stats1(i,j,3,time) = true;
            elseif obs2 == false && mod == false
                % correct negative
                stats1(i,j,4,time) = true;
            end
        end
    end
end

% summmary
hits = sum(stats1(:,:,1,:),4);
miss = sum(stats1(:,:,2,:),4);
fals = sum(stats1(:,:,3,:),4);
corn = sum(stats1(:,:,4,:),4);


% metrics
CSI = hits./(hits + fals + miss);
POD = hits./(hits + miss);
FAR = fals./(hits + fals);



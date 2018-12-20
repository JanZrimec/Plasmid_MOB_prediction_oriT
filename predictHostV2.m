function [host, increp] = predictHostV2(mob,W)
% predictHostV2 Prediction of reportoire of hosts from mob group.

load Host_range_dataV2_30_8_17.mat

if W==1
   host = hostRepMob2{1,mob};
   increp = hostRepMob2{2,mob};
else
   host = hostRepMob200_2{1,mob};
   increp = hostRepMob200_2{2,mob};
end

end
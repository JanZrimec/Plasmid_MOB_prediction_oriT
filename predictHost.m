function [host,increp] = predictHost(mob)
% predictHost Prediction of reportoire of hosts from mob group.

load Host_range_data_19_10_16.mat

host = host_range{mob+1,2};
increp = host_range{mob+1,3};

end
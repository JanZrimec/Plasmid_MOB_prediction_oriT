function [out] = predictHostRange(job,seq,W,T)
% predictHostRange Predict repertoire of hosts of OriT or plasmid.

% predict structural properties - 6 structure model
data = predictStructureOriT(seq);
    
% predict MOB with logistic regression model
[mob, mobnum] = getClass(data,job,W,T);
out.mob = mob; 

% predict host range
[out.host, out.rep] = predictHost(mobnum);
[out.hostV2, out.repV2] = predictHostV2(mobnum,W);

end
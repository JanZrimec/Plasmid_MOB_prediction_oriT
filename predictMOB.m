function [out] = predictMOB(seq,par1,par2)
% wrapper function for predictHostRange.m with additions
% seq ... input oriT sequence 130 bp
% par1 ... model parameter 1 - training data used to construct model
% par2 ... model parameter 2 - amount of structural variables

assert(length(seq)==230, "length(seq) != 230")
assert(par1==64 | par1==200,"par1 != {64,200}")
assert(par2==16 | par2==132,"par2 != {16,132}")
par1_W = [64,200];
par2_T = [132,16];
out = predictHostRange(1,seq,find(par1_W==par1),find(par2_T==par2));

end
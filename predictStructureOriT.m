function [out] = predictStructureOriT(in)
% predictStructureOriT Predict structural properties of OriT region.

tmp = Skripta_napoved_struktur_5_7_2016(in,[10 10 10 10 10]);

% code on interval 10 bp
q=10; 
for i=1:22
    g(i,:) = mean(tmp{1,1}((i-1)*q+1:i*q));
    a(i,:) = mean(tmp{1,2}((i-1)*q+1:i*q));
    h(i,:) = mean(tmp{1,3}((i-1)*q+1:i*q));
    v(i,:) = mean(tmp{1,4}((i-1)*q+1:i*q));
    d(i,:) = mean(tmp{1,5}((i-1)*q+1:i*q));
    t(i,:) = mean(tmp{1,6}((i-1)*q+1:i*q));
end

out{1} = [d' a' h' g' t' v'];

% feature selection
% a = [123;76;80;14;124;98;122;117;120;10;34;79;13;73;78;127];
% out16 = out132(a);

load('ind_16.mat')
out{2} = out{1}(ind_16);

end 
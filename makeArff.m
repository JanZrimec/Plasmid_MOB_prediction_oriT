function [] = makeArff(data,label,name,names)

[m,n] = size(data);
dat2=cell(m,n+1);
%names = ['F','P','Q','V'];

for i = 1:m   
        for j=1:n
            dat2{i,j}=data(i,j);
        end
    dat2{i,end}=names(label(i)); 
end
dlmcell_rand([name,'.txt'],dat2,', ')

fileID = fopen('arff_name_tmp.txt','w');
for i=1:n
fprintf(fileID,'@ATTRIBUTE att%d numeric\n',i);
end
fclose(fileID);

system(['touch arff_name.txt']);
system(['cat arff_name_begin_noname.txt >> arff_name.txt']);
system(['cat arff_name_tmp.txt >> arff_name.txt']);
system(['cat arff_name_end.txt >> arff_name.txt']);
system(['rm arff_name_tmp.txt']);

system(['touch ',name,'.arff']);
system(['cat arff_name.txt >> ',name,'.arff']);
system(['cat ',name,'.txt >> ',name,'.arff']);
system(['rm arff_name.txt']);
system(['rm ',name,'.txt']);

end
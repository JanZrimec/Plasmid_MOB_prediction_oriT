%% predict structure 132 orit 230 bp
% potrebno na novo preracunati set 132 zaradi dolzine in tidd (230bp)

seqs = fastaread('ORIall_204_230bp_sort.fasta');

for i = 1:204
i
    seq = seqs(i).Sequence;
    out = Skripta_napoved_struktur_5_7_2016(seq,[10 10 10 10 10]);
    Pdata_NEW_cell_raw{i,1} = out;
end

q=10; 
for j=1:204
    for i=1:22
        Pdata_NEW_MQ_raw{j,1}(i)=mean(Pdata_NEW_cell_raw{j,1}{1,1}((i-1)*q+1:i*q));   %%% 3.10.2016 - NAPAKA!!! 
        Pdata_NEW_MQ_raw{j,2}(i)=mean(Pdata_NEW_cell_raw{j,1}{1,2}((i-1)*q+1:i*q));
        Pdata_NEW_MQ_raw{j,3}(i)=mean(Pdata_NEW_cell_raw{j,1}{1,3}((i-1)*q+1:i*q));
        Pdata_NEW_MQ_raw{j,4}(i)=mean(Pdata_NEW_cell_raw{j,1}{1,4}((i-1)*q+1:i*q));
        Pdata_NEW_MQ_raw{j,5}(i)=mean(Pdata_NEW_cell_raw{j,1}{1,5}((i-1)*q+1:i*q));
        Pdata_NEW_MQ_raw{j,6}(i)=mean(Pdata_NEW_cell_raw{j,1}{1,6}((i-1)*q+1:i*q));
    end
end

% celoten dataset
for i=1:204

  d(i,:) = Pdata_NEW_MQ_raw{i,5}; %probam z vektorji lastnosti ce lahko
  a(i,:) = Pdata_NEW_MQ_raw{i,2};
  h(i,:) = Pdata_NEW_MQ_raw{i,3};
  g(i,:) = Pdata_NEW_MQ_raw{i,1};
  t(i,:) = Pdata_NEW_MQ_raw{i,6};
  v(i,:) = Pdata_NEW_MQ_raw{i,4};
  
end

Pdata_NEW_raw = [d a h g t v];

[D1, W_lda1] = lda(Pdata_NEW_raw,Pclass_sort);
infoZ1 = cumsum(D1) / sum(D1);
X_proj = project(Pdata_NEW_raw, W_lda1(:,1:3));

function dG = G_NN(sekv) %neki

% 20.1.13

% naloada parametre
load('nn_therm_text.mat') %1:13 (10+3 dod) so wc, 14:61 (48) so sin_int_m, 62:93 (32) so sin_dang_ends
load('nn_therm_data.mat')
load('nn_therm_loops.mat')

%sekv(1,end+1:end+width-1)=sekv(1,1:width-1); % dodatek za cirkularizacijo!!

sekv(2,:)=seqcomplement(sekv);

[~,n] = size(sekv);
NN = [char];
pos = [];

dG = zeros(n-1,1);
%Gw = zeros(n-width+1,1);

for i=1:n-1
    
    NN(i,:) = [sekv(1,i:i+1),'/',sekv(2,i:i+1)];
    pos = strmatch(NN(i,:),nn_therm_text); % poiscemo cez wc, single mismatche, single dangling ende
    
    if isempty(pos) == 1    % ce NN niso pravilno obrnjeni in ne najde
        NN(i,:) = seqreverse(NN(i,:));
        pos = strmatch(NN(i,:),nn_therm_text);
        NN(i,:) = seqreverse(NN(i,:)); %obrne nazaj!
    end  
    
    dG(i) = nn_therm_data(pos,3);
    
end    

% obracun
% for i=1:n-width+1
% 
%     Gw(i) = sum(dG(i:i+width-2)); % zato ker width je stevilo nukleotideov
% 
% end
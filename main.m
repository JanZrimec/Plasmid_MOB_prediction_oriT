%% PredictMOB

function [out] = main(seq,par1,par2)
% wrapper function for predictHostRange.m with additions
% seq ... input oriT sequence 130 bp
% par1 ... model parameter 1 - training data used to construct model
% par2 ... model parameter 2 - amount of structural variables

assert(length(seq)==230, "length(seq) != 230")
%assert(par1==64 | par1==200,"par1 != {64,200}")
%assert(par2==16 | par2==132,"par2 != {16,132}")
par1_W = [64,200];
par2_T = [132,16];
out = predictHostRange(1,seq,find(par1_W==par1),find(par2_T==par2));

end

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

function [mob, mobnum] = getClass(data,job,W,T)
% getClass Predict which MOB group OriT belongs to with pretrained model.

load('model_names.mat')

name1 = model_names{W,T};
name2 = sprintf('job_%d_test',job);
name3 = sprintf('job_%d_out.csv',job);

makeArff(data{T},1,name2,'FPQV')

command = ['java -cp weka.jar weka.classifiers.functions.MultilayerPerceptron -l ',name1,' -T ',name2,'.arff -classifications "weka.classifiers.evaluation.output.prediction.CSV -file ',name3,'"'];
system(command)

tmp = csv2cell(name3,'fromfile');
mobnum = str2num(tmp{2,3}(1));
mob = tmp{2,3}(3);
disp(mob)
disp(mobnum)
system(['rm ',name2,'.arff ',name3])

end

function [host,increp] = predictHost(mob)
% predictHost Prediction of reportoire of hosts from mob group.

load Host_range_data_19_10_16.mat

host = host_range{mob+1,2};
increp = host_range{mob+1,3};

end

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

function [outs]=Skripta_napoved_struktur_5_7_2016(seq,width) %width je vektor 5
% poberem iz nnx predelave za p64

load('NN_structural_properties.mat')
load('NNN_structural_properties.mat')
x = length(seq);
ntseq = nt2int(seq);
jnn = [1 33 34 37];
out = cell(6,1);

for i = 1:4 % specific properties

  j = jnn(i);
  
  for k = 1:x-1  % NN    

      out{i}(k) = NNsp{j+2, 6+(4*(ntseq(k)-1)+ntseq(k+1))};
  end
        
  for k=1:x-(width(i)-1)-1 % sum over width
      if j == 33 %persistence
      
        outs{i}(k) = width(i)/sum(1./out{i}(k:k+width(i)-1));
      
      else
      
        outs{i}(k) = sum(out{i}(k:k+width(i)-1));
        
      end
   end
end
 
for i = 1 % dnaze
    
  for k=1:x-2 % NNN
  
      l=1;
      while (strcmp(NNNsp{l,1}(1:3),seq(k:k+2)) || strcmp(NNNsp{l,1}(5:7),seq(k:k+2))) == 0; %loop breaks when NNN is found
            l=l+1;
      end
      out{5}(k) = NNNsp{l,2};
  end
    
  for k=1:x-(width(5)-2)-2

        outs{5}(k) = sum(out{5}(k:k+width(5)-2));
  end
        
end
    
tmpt = weka_run_klas_10(1,seq,x);
outs{6} = tmpt';

%dG
%persistence
%helical
%deformability
%DNAze
%TIDD

end

function out=weka_run_klas_10(cnn,seq,size)

l=1;
rnd=1;
summ=1;
pbd=3;
wide=6;
tresh=1;
W=10; %:5:20

name = PBDNN_make_arff_var2608_tresh_dud(cnn,W,wide,rnd,summ,pbd,tresh);

bla=sprintf('tmp_%d',1);
name2 = PBDNN_make_arff_euk_nerand(bla,seq,zeros(1,size),W,wide,pbd,tresh);
name3 = [name2,'.test'];
command = ['java -cp weka.jar weka.classifiers.trees.M5P -t ',name,' -T ',name2,' -classifications "weka.classifiers.evaluation.output.prediction.CSV -file ',name3,'"'];
[~]=dos(command);
output=dlmread(name3,',',1,0);
    
out=output(:,3);

end

function name=PBDNN_make_arff_var2608_tresh_dud(cnn,W,wide,rand,summ,p,t)

format long
w=W/5;

switch summ 
case 0

    w1=floor(wide/10);
    w2=wide-w1*10;

    switch p
        case 3
            
            switch rand
            case 0
            name=sprintf('count_W%d_T%d_w%d%d.arff',W,t,w1,w2);
            %dlmcell(name,G,', ')
            case 1
            name=sprintf('count_W%d_T%d_w%d%d_rand.arff',W,t,w1,w2);
            %dlmcell_rand(name,G,', ')
            end
            
        case 4
            
        switch rand
        case 0
            name=sprintf('norm_W%d_T%d_w%d%d.arff',W,t,w1,w2);
            %dlmcell(name,G,', ')
        case 1
            name=sprintf('norm_W%d_T%d_w%d%d_rand.arff',W,t,w1,w2);
            %dlmcell_rand(name,G,', ')
        end
            
        case 5
        
        switch rand
        case 0
            name=sprintf('prob_W%d_T%d_w%d%d.arff',W,t,w1,w2);
            %dlmcell(name,G,', ')
        case 1
            name=sprintf('prob_W%d_T%d_w%d%d_rand.arff',W,t,w1,w2);
            %dlmcell_rand(name,G,', ')
        end
                
    end
    
case 1
   
    w1=floor(wide/10);
    w2=wide-w1*10;

    switch p
        case 3
            
            switch rand
            case 0
            name=sprintf('count_W%d_sum_T%d_w%d%d.arff',W,t,w1,w2);
            %dlmcell(name,G,', ')
            case 1
            name=sprintf('count_W%d_sum_T%d_w%d%d_rand.arff',W,t,w1,w2);
            %dlmcell_rand(name,G,', ')
            end
            
        case 4
            
        switch rand
        case 0
            name=sprintf('norm_W%d_sum_T%d_w%d%d.arff',W,t,w1,w2);
            %dlmcell(name,G,', ')
        case 1
            name=sprintf('norm_W%d_sum_T%d_w%d%d_rand.arff',W,t,w1,w2);
            %dlmcell_rand(name,G,', ')
        end
            
        case 5
        
        switch rand
        case 0
            name=sprintf('prob_W%d_sum_T%d_w%d%d.arff',W,t,w1,w2);
            %dlmcell(name,G,', ')
        case 1
            name=sprintf('prob_W%d_sum_T%d_w%d%d_rand.arff',W,t,w1,w2);
            %dlmcell_rand(name,G,', ')
        end
                
    end

end

end

function name=PBDNN_make_arff_euk_nerand(euk,sekv,input,W,wide,p,tresh)

% 11.3.13 - bomo ubrali cirkularno pot - tako kot pbd rezultati
% 15_3_13 - popravek - napačen preračun

% dodatek za cirkularizacijo!!
n=length(sekv);
sekv=[sekv(1,end-wide+1:end),sekv(1,1:end),sekv(1,1:W-1+wide)];

% rand=1;
% summ=1;
w=W/5;
t=tresh;

% širina wide = wide, medtem ko, ker je bubble vmes, je širina w*5-1 (nn fore pač)

G=cell(n,2+2*wide);


    dG = G_NN(sekv); %je premaknjen in vključuje prej wide in pol wide + width-1
    PBD = input;
        
    for i=1:n
    
        G{i,wide+1}=sum(dG(i-1+wide+1:i-1+wide+W-1));
        
        for l=1:wide
            
            G{i,l}=dG(i-1+l);
            G{i,wide+1+l}=dG(i-1+wide+W-1+l);
        
        end
        
        G{i,2*wide+2}=PBD(i);
        
    end
    
    w1=floor(wide/10);
    w2=wide-w1*10;

    switch p
        
        case 3
            
            name=sprintf('count_%s_W%d_T%d_w%d%d_rand.arff',euk,W,t,w1,w2);
            dlmcell(name,G,', ')
            
        case 4
            
            name=sprintf('norm_%s_W%d_T%d_w%d%d_rand.arff',euk,W,t,w1,w2);
            dlmcell(name,G,', ')
                    
        case 5

            name=sprintf('prob_%s_W%d_T%d_w%d%d_rand.arff',euk,W,t,w1,w2);
            dlmcell(name,G,', ')
                
    end

    text=['@RELATION PBDNN' 10];
    for i=1:wide+1
   
    curr=sprintf('@ATTRIBUTE NN%d	numeric',i-wide);
    text=[text 10 curr];
    
    end
    for i=1:wide
   
    curr=sprintf('@ATTRIBUTE NN%d	numeric',i+W-1);
    text=[text 10 curr];
    
    end
    text=[text 10 '@ATTRIBUTE class	numeric' 10 10 '@DATA'];


dlmwrite(name,[text 13 10 fileread(name)],'delimiter','');

end

function dlmcell(file,cell_array,varargin)
%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %%
% <><><><><>     dlmcell - Write Cell Array to Text File      <><><><><> %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %
%                                                 Version:    01.06.2010 %
%                                                     (c) Roland Pfister %
%                                             roland_pfister@t-online.de %
% 1. Synopsis                                                            %
%                                                                        %
% A single cell array is written to an output file. Cells may consist of %
% any combination of (a) numbers, (b) letters, or (c) words. The inputs  %
% are as follows:                                                        %
%                                                                        %
%       - file       The output filename (string).                       %
%       - cell_array The cell array to be written.                       %
%       - delimiter  Delimiter symbol, e.g. ',' (optional;               %
%                    default: tab ('\t'}).                               %
%       - append     '-a' for appending the content to the               %
%                    output file (optional).                             %
%                                                                        %
% 2. Example                                                             %
%                                                                        %
%         mycell = {'Numbers', 'Letters', 'Words','More Words'; ...      %
%                    1, 'A', 'Apple', {'Apricot'}; ...                   %
%                    2, 'B', 'Banana', {'Blueberry'}; ...                %
%                    3, 'C', 'Cherry', {'Cranberry'}; };                 %
%         dlmcell('mytext.txt',mycell);                                  %
%                                                                        %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %



%% Check input arguments
if nargin < 2
    disp('Error - Give at least two input arguments!'); 
    return;
elseif nargin > 4
    disp('Error - Do not give more than 4 input arguments!'); 
    return;
end
if ~ischar(file)
    disp(['Error - File input has to be a string (e.g. ' ...
        char(39) 'output.txt' char(39) '!']); 
    return;
end;
if ~iscell(cell_array)
    disp('Error - Input cell_array not of the type "cell"!'); 
    return;
end;
delimiter = '\t';
append = 'w';
if nargin > 2
    for i = 1:size(varargin,2)
        if strcmp('-a',varargin{1,i}) == 1
            append = 'a';
        else
            delimiter = varargin{1,i};
        end;
    end;
end

%% Open output file and prepare output array.
output_file = fopen(file,append);
output = cell(size(cell_array,1),size(cell_array,2));

%% Evaluate and write input array.
for i = 1:size(cell_array,1)
for j = 1:size(cell_array,2)
    if numel(cell_array{i,j}) == 0
        output{i,j} = '';
    % Check whether the content of cell i,j is
    % numeric and convert numbers to strings.
    elseif isnumeric(cell_array{i,j}) || islogical(cell_array{i,j})
        output{i,j} = num2str(cell_array{i,j}(1,1));
    
    % Check whether the content of cell i,j is another cell (e.g. a
    % string of length > 1 that was stored as cell. If cell sizes 
    % equal [1,1], convert numbers and char-cells to strings.
    %
    % Note that any other cells-within-the-cell will produce errors
    % or wrong results.
    elseif iscell(cell_array{i,j})
        if size(cell_array{i,j},1) == 1 && size(cell_array{i,j},1) == 1
            if isnumeric(cell_array{i,j}{1,1})
                output{i,j} = num2str(cell_array{i,j}{1,1}(1,1));
            elseif ischar(cell_array{i,j}{1,1})
                 output{i,j} = cell_array{i,j}{1,1};
            end;
        end;
        
     % If the cell already contains a string, nothing has to be done.
     elseif ischar(cell_array{i,j})
         output{i,j} = cell_array{i,j};
     end;
     
     % Cell i,j is written to the output file. A delimiter is appended
     % for all but the last element of each row.
     fprintf(output_file,'%s',output{i,j});
     if j ~= size(cell_array,2)
         fprintf(output_file,'%s',delimiter);
     end
end;
% At the end of a row, a newline is written to the output file.
fprintf(output_file,'\r\n');
end;

%% Close output file.    
fclose(output_file);

end

function dlmcell_rand(file,cell_array,varargin)

precision='%0.9f';

%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %%
% <><><><><>     dlmcell - Write Cell Array to Text File      <><><><><> %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %
%                                                 Version:    01.06.2010 %
%                                                     (c) Roland Pfister %
%                                             roland_pfister@t-online.de %
% 1. Synopsis                                                            %
%                                                                        %
% A single cell array is written to an output file. Cells may consist of %
% any combination of (a) numbers, (b) letters, or (c) words. The inputs  %
% are as follows:                                                        %
%                                                                        %
%       - file       The output filename (string).                       %
%       - cell_array The cell array to be written.                       %
%       - delimiter  Delimiter symbol, e.g. ',' (optional;               %
%                    default: tab ('\t'}).                               %
%       - append     '-a' for appending the content to the               %
%                    output file (optional).                             %
%                                                                        %
% 2. Example                                                             %
%                                                                        %
%         mycell = {'Numbers', 'Letters', 'Words','More Words'; ...      %
%                    1, 'A', 'Apple', {'Apricot'}; ...                   %
%                    2, 'B', 'Banana', {'Blueberry'}; ...                %
%                    3, 'C', 'Cherry', {'Cranberry'}; };                 %
%         dlmcell('mytext.txt',mycell);                                  %
%                                                                        %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %



%% Check input arguments
if nargin < 2
    disp('Error - Give at least two input arguments!'); 
    return;
elseif nargin > 4
    disp('Error - Do not give more than 4 input arguments!'); 
    return;
end
if ~ischar(file)
    disp(['Error - File input has to be a string (e.g. ' ...
        char(39) 'output.txt' char(39) '!']); 
    return;
end;
if ~iscell(cell_array)
    disp('Error - Input cell_array not of the type "cell"!'); 
    return;
end;
delimiter = '\t';
append = 'w';
if nargin > 2
    for i = 1:size(varargin,2)
        if strcmp('-a',varargin{1,i}) == 1
            append = 'a';
        else
            delimiter = varargin{1,i};
        end;
    end;
end

%% Open output file and prepare output array.
output_file = fopen(file,append);
output = cell(size(cell_array,1),size(cell_array,2));

%% Evaluate and write input array.
for i = randperm(size(cell_array,1))
for j = 1:size(cell_array,2)
    if numel(cell_array{i,j}) == 0
        output{i,j} = '';
    % Check whether the content of cell i,j is
    % numeric and convert numbers to strings.
    elseif isnumeric(cell_array{i,j}) || islogical(cell_array{i,j})
        output{i,j} = num2str(cell_array{i,j}(1,1),precision);
    
    % Check whether the content of cell i,j is another cell (e.g. a
    % string of length > 1 that was stored as cell. If cell sizes 
    % equal [1,1], convert numbers and char-cells to strings.
    %
    % Note that any other cells-within-the-cell will produce errors
    % or wrong results.
    elseif iscell(cell_array{i,j})
        if size(cell_array{i,j},1) == 1 && size(cell_array{i,j},1) == 1
            if isnumeric(cell_array{i,j}{1,1})
                output{i,j} = num2str(cell_array{i,j}{1,1}(1,1),precision);
            elseif ischar(cell_array{i,j}{1,1})
                 output{i,j} = cell_array{i,j}{1,1};
            end;
        end;
        
     % If the cell already contains a string, nothing has to be done.
     elseif ischar(cell_array{i,j})
         output{i,j} = cell_array{i,j};
     end;
     
     % Cell i,j is written to the output file. A delimiter is appended
     % for all but the last element of each row.
     fprintf(output_file,'%s',output{i,j});
     if j ~= size(cell_array,2)
         fprintf(output_file,'%s',delimiter);
     end
end;
% At the end of a row, a newline is written to the output file.
fprintf(output_file,'\r\n');
end;

%% Close output file.    
fclose(output_file);

end

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
end

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

function data = csv2cell(varargin)
% CSV2CELL - parses a Windows CSV file into an NxM cell array, where N is
% the number of lines in the CSV text and M is the number of fields in the
% longest line of the CSV file. Lines are delimited by carriage returns
% and/or newlines.
%
% A Windows CSV file format allows for commas (,) and double quotes (") to
% be contained within fields of the CSV file. Regular fields are just text
% separated by commas (e.g. foo,bar,hello world). Fields containing commas
% or double quotes are surrounded by double quotes (e.g.
% foo,bar,"item1,item2,item3",hello world). In the previous example,
% "item1,item2,item3" is one field in the CSV text. For double quotes to be
% represented, they are written in pairs in the file, and contained within
% a quoted field, (e.g. foo,"this field contains ""quotes""",bar). Spaces
% within fields (even leading and trailing) are preserved.
%
% All fields from the CSV file are returned as strings. If the CSV text
% contains lines with different numbers of fields, then the "missing"
% fields with appear as empty arrays, [], in the returned data. You can
% easily convert the data you expect to be numeric using str2num() and
% num2cell(). 
%
% Examples:
%  >> csv2cell('foo.csv','fromfile') % loads and parses entire file
%  >> csv2cell(',,,') % returns cell array {'','','',''}
%  >> csv2cell(',,,','text') % same as above, declaring text input
%  >> csv2cell(sprintf('%s\r\n%s',...
%     '"Ten Thousand",10000,,"10,000","""It''s ""10 Grand"", baby",10k',...
%     ',foo,bar,soo'))
%  ans = 
%    'Ten Thousand'    '10000'       ''    '10,000'    [1x22 char]    '10k'
%                ''    'foo'      'bar'    'soo'                []       []
%  >> % note the two empty [] cells, because the second line has two fewer
%  >> % fields than the first. The empty field '' at the beginning of the
%  >> % second line is due to the leading comma on that line, which is
%  >> % correct behavior. A trailing comma will do the same to the end of a
%  >> % line.
% 
% Limitations/Exceptions:
%   * This code is untested on large files. It may take a long time due to
%   variables growing inside loops (yes, poor practice, but easy coding).
%   * This code has been minimally tested to work with a variety of weird
%   Excel files that I have.
%   * Behavior with improperly formatted CSV files is untested.
%   * Technically, CSV files from Excel always separate lines with the pair
%   of characters \r\n. This parser will also separate lines that have only
%   \r or \n as line terminators. 
%   * Line separation is the first operation. I don't think the Excel CSV
%   format has any allowance for newlines or carriage returns within
%   fields. If it does, then this parser does not support it and would not
%   return bad output.
%
% Copyright 2008 Arthur Hebert

% Process arguments
if nargin == 1
    text = varargin{1};
elseif nargin == 2
    switch varargin{2}
        case 'fromfile'
            filename = varargin{1};
            fid = fopen(filename);
            text = char(fread(fid))';
            fclose(fid);
        case 'text'
            text = varargin{1};
        otherwise
            error('Invalid 2nd argument %s. Valid options are ''fromfile'' and ''text''',varargin{2})
    end
else
    error('CSV2CELL requires 1 or 2 arguments.')
end


% First split it into lines
lines = regexp(text,'(\r\n|[\r\n])','split'); % lines should now be a cell array of text split by newlines

% a character is either a delimiter or a field
inField = true;
inQuoteField = false;
% if inField && ~inQuoteField --> then we're in a raw field

skipNext = false;
data = {};
field = '';
for lineNumber = 1:length(lines)
    nChars = length(lines{lineNumber}); % number of characters in this line
    fieldNumber = 1;
    for charNumber = 1:nChars
        if skipNext
            skipNext = false;
            continue
        end
        thisChar = lines{lineNumber}(charNumber);
        if thisChar == ','
            if inField 
                if inQuoteField % this comma is part of the field
                    field(end+1) = thisChar;
                else % this comma is the delimiter marking the end of the field
                    data{lineNumber,fieldNumber} = field;
                    field = '';
                    fieldNumber = fieldNumber + 1;
                end
            else % we are not currently in a field -- this is the start of a new delimiter
                inField = true;
            end
            if charNumber == nChars % this is a hanging comma, indicating the last field is blank
                data{lineNumber,fieldNumber} = '';
                field = '';
                fieldNumber = fieldNumber + 1;
            end
        elseif thisChar == '"' 
            if inField
                if inQuoteField
                    if charNumber == nChars % it's the last character, so this must be the closing delimiter?
                        inField = false;
                        inQuoteField = false;
                        data{lineNumber,fieldNumber} = field;
                        field = '';
                        fieldNumber = fieldNumber + 1;
                    else 
                        if lines{lineNumber}(charNumber+1) == '"' % this is translated to be a double quote in the field
                            field(end+1) = '"';
                            skipNext = true;
                        else % this " is the delimiter ending this field
                            data{lineNumber,fieldNumber} = field;
                            field = '';
                            inField = false;
                            inQuoteField = false;
                            fieldNumber = fieldNumber + 1;
                        end
                    end
                else % this is a delimiter and we are in a new quote field
                    inQuoteField = true;
                end
            else % we are not in a field. This must be an opening quote for the first field?
                inField = true;
                inQuoteField = true;
            end
        else % any other character ought to be added to field
            field(end+1) = thisChar;
            if charNumber == nChars
                data{lineNumber,fieldNumber} = field;
                field = '';
                fieldNumber = fieldNumber + 1;
            elseif charNumber == 1 % we are starting a new raw field
                inField = true;
            end
        end
    end
end
            
end

function [out] = seqcomplement(in)

in = upper(in);
n = length(in);

for i = 1:n
  switch in(i)
  case 'A'
    out(i) = 'T';
  case 'C'
    out(i) = 'G';
  case 'G'
    out(i) = 'C';
  case 'T'
    out(i) = 'A';
  otherwise 
    out(i) = 'A'; %foolproof
  end
  
end
end

function [out] = seqreverse(in)

n = length(in);
out=in(n:-1:1)

end

function [out] = nt2int(in)

in = upper(in);
n = length(in);

for i = 1:n
  switch in(i)
  case 'A'
    out(i) = 1;
  case 'C'
    out(i) = 2;
  case 'G'
    out(i) = 3;
  case 'T'
    out(i) = 4;
  otherwise 
    out(i) = 5;
  end
  
end
end

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


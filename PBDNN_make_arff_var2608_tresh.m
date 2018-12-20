function name=PBDNN_make_arff_var2608_tresh(cnn,W,wide,rand,summ,p,t)

% 2.3.13 - all in one
% all variable - sum = 2608
format long
w=W/5;

% širina wide = wide, medtem ko, ker je bubble vmes, je širina w*5-1 (nn fore pač)

no=10; %št.plazmidov
st=30; %

switch summ
    
case 0

    G=cell(4*340+315+4*351+329-10*(2*st+20),W-1+2*wide+1); %2608 %PBD norm, P
    k=0;

    for i=1:10

    sekv = cnn{2,i+1};
    dG = G_NN(sekv);
    PBD = cnn{p,i+1}(:,w,t);
    [ed,~,~] = size(cnn{p,i+1});
    ed = ed-st-20;
    
    for j=1:ed-st %začetek zamaknjen, tako kot konec, za upoštevanje wide %235
        
        for l=1:W-1+2*wide
            
            G{k+j,l}=dG(j+l-1);
        
        end    
        
        G{k+j,W-1+2*wide+1}=PBD(j+wide);
        
    end
    
    k=k+ed-st;
    end

    w1=floor(wide/10);
    w2=wide-w1*10;

    switch p
        case 3
            
            switch rand
            case 0
            name=sprintf('count_W%d_T%d_w%d%d.arff',W,t,w1,w2);
            dlmcell(name,G,', ')
            case 1
            name=sprintf('count_W%d_T%d_w%d%d_rand.arff',W,t,w1,w2);
            dlmcell_rand(name,G,', ')
            end
            
        case 4
            
        switch rand
        case 0
            name=sprintf('norm_W%d_T%d_w%d%d.arff',W,t,w1,w2);
            dlmcell(name,G,', ')
        case 1
            name=sprintf('norm_W%d_T%d_w%d%d_rand.arff',W,t,w1,w2);
            dlmcell_rand(name,G,', ')
        end
            
        case 5
        
        switch rand
        case 0
            name=sprintf('prob_W%d_T%d_w%d%d.arff',W,t,w1,w2);
            dlmcell(name,G,', ')
        case 1
            name=sprintf('prob_W%d_T%d_w%d%d_rand.arff',W,t,w1,w2);
            dlmcell_rand(name,G,', ')
        end
                
    end

    text=['@RELATION PBDNN' 10];
    for i=1:W-1+2*wide
   
    curr=sprintf('@ATTRIBUTE NN%d	numeric',i-wide);
    text=[text 10 curr];
    
    end
    text=[text 10 '@ATTRIBUTE class	numeric' 10 10 '@DATA'];

case 1
        
    G=cell(4*340+315+4*351+329-10*(2*st+20),2*wide+2); %2608 %PBD norm, P
    k=0;

    for i=1:10

    sekv = cnn{2,i+1};
    dG = G_NN(sekv);
    PBD = cnn{p,i+1}(:,w,t);
    [ed,~,~] = size(cnn{p,i+1});
    ed = ed-st-20;
    
    for j=1:ed-st %začetek zamaknjen, tako kot konec, za upoštevanje wide

        G{k+j,wide+1}=sum(dG(j-1+wide+1:j-1+wide+W-1));
        
        for l=1:wide
            
            G{k+j,l}=dG(j-1+l);
            G{k+j,wide+1+l}=dG(j-1+wide+W-1+l);
        
        end
        
        G{k+j,2*wide+2}=PBD(j+wide);
        
    end
    
    k=k+ed-st;
    end

    w1=floor(wide/10);
    w2=wide-w1*10;

    switch p
        case 3
            
            switch rand
            case 0
            name=sprintf('count_W%d_sum_T%d_w%d%d.arff',W,t,w1,w2);
            dlmcell(name,G,', ')
            case 1
            name=sprintf('count_W%d_sum_T%d_w%d%d_rand.arff',W,t,w1,w2);
            dlmcell_rand(name,G,', ')
            end
            
        case 4
            
        switch rand
        case 0
            name=sprintf('norm_W%d_sum_T%d_w%d%d.arff',W,t,w1,w2);
            dlmcell(name,G,', ')
        case 1
            name=sprintf('norm_W%d_sum_T%d_w%d%d_rand.arff',W,t,w1,w2);
            dlmcell_rand(name,G,', ')
        end
            
        case 5
        
        switch rand
        case 0
            name=sprintf('prob_W%d_sum_T%d_w%d%d.arff',W,t,w1,w2);
            dlmcell(name,G,', ')
        case 1
            name=sprintf('prob_W%d_sum_T%d_w%d%d_rand.arff',W,t,w1,w2);
            dlmcell_rand(name,G,', ')
        end
                
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

end

dlmwrite(name,[text 13 10 fileread(name)],'delimiter','');

end
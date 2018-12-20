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
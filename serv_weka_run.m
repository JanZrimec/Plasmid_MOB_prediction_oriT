function out=serv_weka_run(job,sekv,W,T,wide)

% Main file for prediction of TIDD based on Weka models built by M5P algorithm (see Zrimec, 2014) 

% job	.. job number (for batch predictions)
% sekv	.. nucleotide sequence (capital letters, ATGC only)
% W	.. destabilization width
% T	.. destabilization treshold
% wide	.. number of included neighboring regions (see Figure 1)

n=length(sekv);
name = serv_name(W,T,wide);
name2 = serv_name2(job,sekv,n,W,T,wide);
name3 = sprintf('job%d.out',job);

command = ['java -cp weka.jar weka.classifiers.trees.M5P -t ',name,' -T ',name2,' -classifications "weka.classifiers.evaluation.output.prediction.CSV -file ',name3,'"'];
[~]=system(command);
result=dlmread(name3,',',1,0);
out=result(:,3);

end

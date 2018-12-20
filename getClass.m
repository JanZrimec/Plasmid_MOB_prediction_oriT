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
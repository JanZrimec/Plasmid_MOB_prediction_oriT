makeArff(Pdata_NEW3_raw,y,'data_204_132',names)
makeArff(Pdata_NEW3_raw(:,ind_16),y,'data_204_16',names)
makeArff(Pdata_NEW3_raw(Pindex_64,:),y(Pindex_64),'data_64_132',names)
makeArff(Pdata_NEW3_raw(Pindex_64,ind_16),y(Pindex_64),'data_64_16',names)

system('java -cp weka.jar weka.classifiers.functions.MultilayerPerceptron -t data_204_132.arff -d data_204_132.model');
system('java -cp weka.jar weka.classifiers.functions.MultilayerPerceptron -t data_204_16.arff -d data_204_16.model');
system('java -cp weka.jar weka.classifiers.functions.MultilayerPerceptron -t data_64_132.arff -d data_64_132.model');
system('java -cp weka.jar weka.classifiers.functions.MultilayerPerceptron -t data_64_16.arff -d data_64_16.model');

seqs = fastaread('ORIall_204_230bp_sort.fasta');
in = seq(204).Sequence;
out = predictStructureOriT(in);

makeArff(out,1,'test_seq',names)

name1 = 'data_204_132.model';
name2 = 'job_1_test.arff';
name3 = 'job_1_out.csv';

command = ['java -cp weka.jar weka.classifiers.functions.MultilayerPerceptron -l ',name1,' -T ',name2,' -classifications "weka.classifiers.evaluation.output.prediction.CSV -file ',name3,'"'];
system(command)

tmp = csv2cell(name3,'fromfile');
mobnum = tmp{2,3}(1);
mob = tmp{2,3}(3);

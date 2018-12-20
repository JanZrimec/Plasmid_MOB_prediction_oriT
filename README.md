# Prediction of plasmid MOB groups based on DNA structures in *oriT*

Tools and data from [Zrimec & Lapanje 2018: DNA structure at the plasmid origin-of-transfer indicates its potential transfer range](https://www.nature.com/articles/s41598-018-20157-y)

To reproduce results from the paper additional source code is available upon request.

## Description
Horizontal gene transfer via plasmid conjugation enables antimicrobial resistance (AMR) to spread among bacteria and is a major health concern. 

Increased incidence of antimicrobial resistance (AMR) in bacteria has raised global awareness in the recent years and will call for even more attention in the following decades, due to the decreased rate of introduction of new antibiotics. One of the approaches to treat infections of AMR bacteria is to use the most appropriate antibiotic combinations and treatment regimens that inhibit the transfer of AMR genes, since inappropriate treatments increase the number of AMR variants. In human pathogenic bacteria, AMR genes have arisen through transfer of mobile elements from the large AMR environmental pool. One of the most important of these elements are conjugative plasmids, where each plasmid can be transferred to and hosted in a particular repertoire of appropriate host bacteria. This is characterized by the plasmids MOB group. Given that the potential host range of a plasmid is not defined only by plasmid transfer, but also by the propensity of the plasmid to stabilize in the subsequent generations of the bacterial host, two separate ranges can be distinguished: (i) the range of potential transfer hosts, based on the hosts of plasmids used for training the models (see paper Supp. Table S7), and (ii) the range of potential incompatibility and replication (Inc/Rep) types that can help determine the replication host range (see paper Supp. Table S8).

<img src="https://github.com/JanZrimec/Plasmid_MOB_prediction_oriT/blob/master/Figure1.png" width="480">

To facilitate prediction of plasmid MOB groups, we have developed a bioinformatic procedure based on analysis of the origin-of-transfer (*oriT*), a merely 230 bp long non-coding plasmid DNA region that is the enzymatic substrate for the relaxase. By computationally interpreting conformational and physicochemical properties of the oriT region, which facilitate relaxase-*oriT* recognition and initiation of nicking, MOB groups can be resolved with over 99% accuracy. We have shown that *oriT* structural properties are highly conserved and can be used to discriminate among MOB groups more efficiently than the *oriT* nucleotide sequence. The procedure for prediction of MOB groups and potential transfer range of plasmids can also be run at http://dnatools.eu.

## Usage
```out = predictMOB(seq,par1,par2)```

where:
* seq ... input *oriT* sequence of length 230 bp type char
* par1 ... model parameter 1: training data used to construct model: {64, 200}
* par2 ... model parameter 2: amount of structural variables: {16, 132}

* out ... structure array with results:
  * out.mob ... predicted MOB group
  * out.host ... refer to out.hostV2 for results published in final revision of paper
  * out.rep ... refer to out.repV2 for results published in final revision of paper
  * out.hostV2 ... ranges of potential transfer hosts based on pooling known transfer host clades found in the training sets of either 64 or 200 elements
  * out.repV2 ... ranges of potential incopatibility and replication (Inc/Rep) types based on pooling known Inc/Rep types found in the training sets of either 64 or 200 elements.

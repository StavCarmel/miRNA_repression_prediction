clear all;
close all;
clc

%% What should this script do
create_training_set=0;  %turn on (1) for finding sites in the training set
create_model=0;         %turn on (1) for building a regressor
create_validation_set=0;  %turn on (1) for finding sites in the validation set
predict_val=1;          %turn on (1) for predicting the val dataset

%% Reading relevant data
load('genes_training.mat')
load('miRs_training.mat') 
load('repress.mat')
miRs_seq=values(miRs);
miRs_names=keys(miRs);

%% Finds sites in training set
if create_training_set
    [training_ORF,training_UTR3] = extract_bindsites(genes,miRs,repress,1);
    save('training_ORF','training_ORF');
    save('training_UTR3','training_UTR3');
else
    load('training_ORF.mat');
    load('training_UTR3.mat');
end

%% Creates model
if create_model
    for i=1:height(training_ORF)
        genes_num=training_ORF.gene(i);
        feat_ORF_train(i,:)=extract_features(training_ORF(i,:),genes(genes_num,:),'ORF');
    end

    for j=1:height(training_UTR3)
        genes_num=training_UTR3.gene(j);
        feat_UTR3_train(j,:)=extract_features(training_UTR3(j,:),genes(genes_num,:),'UTR3');
    end
    
    [B,FitInfo]=lasso(feat_ORF_train,training_ORF.repress,'CV',10,'Alpha',0.5);
    indORF=FitInfo.IndexMinMSE; %chooses lambdas with min MSE to training set
    ORF_model=[FitInfo.Intercept(indORF);B(:,indORF)];
    
    [B,FitInfo]=lasso(feat_UTR3_train,training_UTR3.repress,'CV',10,'Alpha',0.5);
    indUTR3=FitInfo.IndexMinMSE; %chooses lambdas with min MSE to training set
    UTR3_model=[FitInfo.Intercept(indUTR3);B(:,indUTR3)];

    save('regressor_models','ORF_model','UTR3_model');
    save('training_features','feat_ORF_train','feat_UTR3_train');
else
    load('regressor_models');
end

%% Finds sites in validation set
load('genes_validation.mat')
load('miR_validation.mat') 
if create_validation_set
    [validation_ORF,validation_UTR3] = extract_bindsites(genes,miRs,0,0);
    save('validation_ORF','validation_ORF');
    save('validation_UTR3','validation_UTR3');
else
    load('validation_ORF');
    load('validation_UTR3');
end

%% Predicts validation's repression
if predict_val
    for k=1:height(validation_ORF)
        genes_num=validation_ORF.gene(k);
        feat_ORF_val(k,:)=extract_features(validation_ORF(k,:),genes(genes_num,:),'ORF');
    end
    feat_ORF_val=[ones(length(feat_ORF_val),1),feat_ORF_val];

    for m=1:height(validation_UTR3)
        genes_num=validation_UTR3.gene(m);
        feat_UTR3_val(m,:)=extract_features(validation_UTR3(m,:),genes(genes_num,:),'UTR3');
    end
    feat_UTR3_val=[ones(length(feat_UTR3_val),1),feat_UTR3_val];

    predicted_ORF=sum(feat_ORF_val.*ORF_model',2);
    predicted_UTR3=sum(feat_UTR3_val.*UTR3_model',2);

    load('miR_validation.mat')
    miRs_names=keys(miRs);

    predicted_per_pair=cell(height(genes),length(miRs_names));

    for gene=1:height(genes)
        predicted_per_pair{gene,1}=genes.ID(gene);
    end

    for gene=1:height(genes)
        for mi=1:length(miRs_names)
            idx_pair_ORF=find(strcmp(validation_ORF.miRNA_name,miRs_names(mi)).*(validation_ORF.gene==gene));
            idx_pair_UTR3=find(strcmp(validation_UTR3.miRNA_name,miRs_names(mi)).*(validation_UTR3.gene==gene));
        
            cur_sum=sum(predicted_ORF(idx_pair_ORF))+sum(predicted_UTR3(idx_pair_UTR3));
            predicted_per_pair{gene,mi+1}=cur_sum;
        end
    end
    predicted_repression_table=cell2table(predicted_per_pair,'VariableNames',{'ID',miRs_names{:}});
    save('predicted_repression_table','predicted_repression_table');
end



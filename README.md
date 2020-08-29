# miRNA_repression_prediction
In this project I predicted gene's repressions due to miRNA (micro-RNA) binding 

Functions:
Challenge_main - the main script that calls all the functions and creates the model.
Extract_bindsites - this function finds the bindsited and creates a training table according to the position on the mRNA (ORF/3'UTR) and according to the dataset (training/validation).
Extract_features - creates the features.

Veriables:
Training_ORF/UTR3 - tables that contains the details on all the bindsites at the ORF/3'UTR, including gene's number, bindsite's index, miRNA name, repression value.
Training_features - includes 'feat_OR_train' and 'feat_UTR3_train' - features matrices for every bindsite, each row is a bindsite and every column is a different feature.
Validation_ORF/UTR3 - tables similar to Training_ORF/UTR3 but without the repression value.
Regressor_models - contains 'UTR3_model' and 'ORE_model' - a vector contains the wheights for every feature (the first value corresponds to ones vector).
Predicted_repression_table - a table of the repression predictions for every gene.

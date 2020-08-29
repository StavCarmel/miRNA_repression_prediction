function feature_vec = extract_features(site_info,gene_info,region)
% string-region, is the binding site region- ORF or UTR3

am_f=15;                             % amount of features
feature_vec=zeros(1,am_f);           % every column is different feature
bind_ind=site_info.index;            % the binding index in the whole gene
miRNA_l=8;                           % length of site
gene_seq=strcat(gene_info.UTR5{1},gene_info.ORF{1},gene_info.UTR3{1});  % gene sequence

%% folding energy with vienna
%----- RNAcofold -----
fid=fopen('seqs.txt','wt'); %creating the txt file
fprintf(fid,[gene_seq(bind_ind:bind_ind+7),'&',site_info.mi_RNA_seq{1}(1:8),'\n@']);
fclose(fid);
[~,output]=system('RNAcofold.exe<seqs.txt'); 
hibrid_score_cofold=(regexp(output,'-?\d+\.\d*','match')); %searching for the energy number

%----- RNAfold ------
% Calculates fold on 20 nts upstream and downstream a 35 flank around the
% site
conditions(1:length(gene_seq))='.'; %creating a dot string
if bind_ind<=35
    left=1;
else
    left=bind_ind-35;
end
if length(gene_seq)<bind_ind+8+35
    right=length(gene_seq);
else
    right=bind_ind+8+35;
end

% right side and left are long enough
if bind_ind+35+8+20<=length(gene_seq) && bind_ind-20-35>0
    conditions=conditions(bind_ind-20-35:bind_ind+35+8+20);
    new_gene_seq=gene_seq(bind_ind-20-35:bind_ind+35+8+20);
% right side isn't long enough but left is long enough
elseif bind_ind+35+8+20>length(gene_seq) && bind_ind-20-35>0
       conditions=conditions(bind_ind-20-35:end);
       new_gene_seq=gene_seq(bind_ind-20-35:end);
% right side is long enough but left isn't long enough
elseif bind_ind-20-35<=0 && bind_ind+35+8+20<=length(gene_seq)
        conditions=conditions(1:bind_ind+35+8+20);
        new_gene_seq=gene_seq(1:bind_ind+35+8+20);
% right and left arn't long enough
elseif bind_ind-20-35<=0 && bind_ind+35+8+20>length(gene_seq)
        new_gene_seq=gene_seq;
end
         
conditions(left:right)='x'; % inserting 'X' in the 8 binding positions in the gene
fid=fopen('seqs.txt','wt');
fprintf(fid,[new_gene_seq,'\n',conditions,'\n@']); 
fclose(fid);
[~,output2]=system('RNAfold.exe < seqs.txt');
hibrid_score_fold=(regexp(output2,'-?\d+\.\d*','match')); %searching for the energy number

feature_vec(1,1)=str2num(hibrid_score_cofold{1});
feature_vec(1,2)=str2num(hibrid_score_fold{1});

%% conservations
conservations=gene_info.conservation{1};

cons_binding_site=mean(conservations(bind_ind:bind_ind+miRNA_l-1)); %mean on the al the binding positins conservations in this specific binding site
feature_vec(1,3)=cons_binding_site;

%% sequences lengths
ORF_l=length(gene_info.ORF{1});
UTR3_l=length(gene_info.UTR3{1});
UTR5_l=length(gene_info.UTR5{1});

feature_vec(1,4)=ORF_l;
feature_vec(1,5)=UTR3_l;
feature_vec(1,6)=UTR5_l;

%% GC content
GC_ORF=(count(gene_info.ORF{1},'G')+count(gene_info.ORF{1},'C'))/length(gene_info.ORF{1});
GC_UTR3=(count(gene_info.UTR3{1},'G')+count(gene_info.UTR3{1},'C'))/length(gene_info.UTR3{1});
GC_UTR5=(count(gene_info.UTR5{1},'G')+count(gene_info.UTR5{1},'C'))/length(gene_info.UTR5{1});

feature_vec(1,7)=GC_ORF;
feature_vec(1,8)=GC_UTR3;
feature_vec(1,9)=GC_UTR5;

%% 7mer-m8 
base_pairs=seqrcomplement(site_info.mi_RNA_seq{1}(2:7));
pos_8=seqrcomplement(site_info.mi_RNA_seq{1}(8));
mer7_ORF=regexp(gene_info.ORF{1},[pos_8,'(?=',base_pairs,'[^A])']);
amount_7mer_ORF=length(mer7_ORF);

mer7_UTR3=regexp(gene_info.UTR3{1},[pos_8,'(?=',base_pairs,'[^A])']);
amount_7mer_UTR3=length(mer7_UTR3);

feature_vec(1,10)=amount_7mer_ORF;
feature_vec(1,11)=amount_7mer_UTR3;

%% distance from edge
%checking what class of binding sites in order to mach the index to the site
if strcmp(region,'ORF')
    bind_ind_ORF=bind_ind-length(gene_info.UTR5{1});
    edge_dist=min([bind_ind_ORF,length(gene_info.ORF{1})-(bind_ind_ORF+miRNA_l-1)]); % +8?
elseif strcmp(region,'UTR3')
       bind_ind_UTR3=bind_ind-length(gene_info.UTR5{1})-length(gene_info.ORF{1});
       edge_dist=min([bind_ind_UTR3,length(gene_info.UTR3{1})-(bind_ind_UTR3+miRNA_l-1)]); %+8?
end
feature_vec(1,12)=edge_dist;

%% CAI & tAI
load('codon_weights.mat');
CAI=[];
tAI=[];
for i=3:length(gene_info.ORF{1})/3
    curr_codon=gene_info.ORF{1}(i*3-2:i*3);
    CAI=[CAI,CAI_weights(curr_codon)];
    tAI=[tAI,tAI_weights(curr_codon)];
end
feature_vec(1,13)=geomean(CAI);
feature_vec(1,14)=geomean(tAI);

%% amino acid charge
aa_ORF=nt2aa(gene_info.ORF{1}); 
positive_aa=count(aa_ORF,{'R','K','H'});
negative_aa=count(aa_ORF,{'E','D'});
aa_cahrge=(positive_aa-negative_aa)/length(gene_info.ORF{1}); %normalization to length
feature_vec(1,15)=aa_cahrge;

end


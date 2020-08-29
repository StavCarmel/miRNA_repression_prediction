function [training_ORF,training_UTR3] = extract_bindsites(genes,miRNAs,repress,istrain)

miRs_seq=values(miRNAs);
miRs_names=keys(miRNAs);

if istrain==1
    % creating a transition matrix for the repression columns: the first
    % miRNA in the map is in the trans_miRs(1) column of repress
    col_repress=repress.Properties.VariableNames;
    repress_cell=table2cell(repress);
    for mi=1:length(miRs_names)
        trans_miRs(mi)=find(strcmp(col_repress(1,:),miRs_names{mi}));
    end
end

match_ind={};
count_ORF=0;
count_UTR3=0;
for gene=1:length(genes.ORF) %runs on genes
    for mi=1:length(miRs_names) %runs on miRNA
        cur_miRs=miRs_seq(mi); %the miRNA sequence
        base_pairs=seqrcomplement(cur_miRs{1}(2:8));
        new_bp = strrep(base_pairs,'U','T'); %change 'U' to 'T'
        mer8=strcat(new_bp,'A');
        index_ORF=strfind(genes.ORF{gene},mer8);
        index_UTR3=strfind(genes.UTR3{gene},mer8);
        orfutr=0;
        if istrain==1
            if (length(index_ORF) + length(index_UTR3))==1 %Only one site found
                if length(index_ORF)==1
                   orfutr=1;
                elseif length(index_UTR3)==1
                   orfutr=2;
                end
                if orfutr==1
                   count_ORF=count_ORF+1;
                   match_ORF{count_ORF,1}=gene;
                   match_ORF{count_ORF,2}=index_ORF+length(genes.UTR5{gene}); % the index in the whole gene
                   match_ORF{count_ORF,3}=miRs_names{mi};
                   match_ORF{count_ORF,4}=miRs_seq{mi};
                   match_ORF{count_ORF,5}=repress_cell{gene,trans_miRs(mi)};
                elseif orfutr==2
                    count_UTR3=count_UTR3+1;
                    match_UTR3{count_UTR3,1}=gene;
                    match_UTR3{count_UTR3,2}=index_UTR3+length(genes.UTR5{gene})+length(genes.ORF{gene});
                    match_UTR3{count_UTR3,3}=miRs_names{mi};
                    match_UTR3{count_UTR3,4}=miRs_seq{mi};
                    match_UTR3{count_UTR3,5}=repress_cell{gene,trans_miRs(mi)};
                end
            end
        else %not the training set
            for i=1:length(index_ORF)
                count_ORF=count_ORF+1;
                match_ORF{count_ORF,1}=gene;
                match_ORF{count_ORF,2}=index_ORF(i)+length(genes.UTR5{gene}); % the index in the whole gene
                match_ORF{count_ORF,3}=miRs_names{mi};
                match_ORF{count_ORF,4}=miRs_seq{mi};
            end
            for i=1:length(index_UTR3)
                count_UTR3=count_UTR3+1;
                match_UTR3{count_UTR3,1}=gene;
                match_UTR3{count_UTR3,2}=index_UTR3(i)+length(genes.UTR5{gene})+length(genes.ORF{gene});
                match_UTR3{count_UTR3,3}=miRs_names{mi};
                match_UTR3{count_UTR3,4}=miRs_seq{mi};
            end
        end
    end
end
if istrain==1
    training_ORF=cell2table(match_ORF,'VariableNames',{'gene','index','miRNA_name','mi_RNA_seq','repress'});
    training_UTR3=cell2table(match_UTR3,'VariableNames',{'gene','index','miRNA_name','mi_RNA_seq','repress'});
else
    training_ORF=cell2table(match_ORF,'VariableNames',{'gene','index','miRNA_name','mi_RNA_seq'});
    training_UTR3=cell2table(match_UTR3,'VariableNames',{'gene','index','miRNA_name','mi_RNA_seq'});
end
end


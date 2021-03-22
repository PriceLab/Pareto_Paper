clear
load('Parti_vars.mat');
load('clinical_labs_ls.mat');

n_arc=4;

corr_coef=zeros(length(analytes), n_arc);
P_coef=zeros(length(analytes), n_arc);

%[coefs1,scores1,variances] = pca(DataPoints);

for i=1:length(analytes)
    for j=1:n_arc
       [tmp, p]=corrcoef(pc(:,j), DataPoints(:,i));
       corr_coef(i,j)=tmp(1,2);
       P_coef(i,j)=p(1,2);
    end
end

[sort_pc1, pc1_idx] = sort(corr_coef(:,1), 'descend');
[sort_pc2, pc2_idx] = sort(corr_coef(:,2), 'descend');
[sort_pc3, pc3_idx] = sort(corr_coef(:,3), 'descend');
[sort_pc4, pc4_idx] = sort(corr_coef(:,4), 'descend');

names_1=analytes(pc1_idx)';
names_2=analytes(pc2_idx)';
names_3=analytes(pc3_idx)';
names_4=analytes(pc4_idx)';

P_1=P_coef(:,1); P_1=P_1(pc1_idx);
P_2=P_coef(:,2); P_2=P_2(pc2_idx);
P_3=P_coef(:,3); P_3=P_3(pc3_idx);
P_4=P_coef(:,4); P_4=P_4(pc4_idx);

all_corr=[sort_pc1, sort_pc2, sort_pc3, sort_pc4];
all_names=[names_1, names_2, names_3, names_4];


fullTable=cell(length(sort_pc1)+1, 8);
Titles={'PC1', 'R_1', 'PC2', 'R_2','PC3','R_3','PC4', 'R_4',};
fullTable(1,:) = Titles;
fullTable(2:end,1) = names_1;
fullTable(2:end,2) = num2cell(sort_pc1);

fullTable(2:end,3) = names_2;
fullTable(2:end,4) = num2cell(sort_pc2);


fullTable(2:end,5) = names_3;
fullTable(2:end,6) = num2cell(sort_pc3);

fullTable(2:end,7) = names_4;
fullTable(2:end,8) = num2cell(sort_pc4);


cell2csv('corr_pc_analyte.csv', fullTable);

clear

load('Parti_vars.mat');
load('clinical_labs_ls.mat');

n_arc=4;

corr_coef=zeros(length(analytes), n_arc);
P_coef=zeros(length(analytes), n_arc);

[~, dists] = sortDataByDistance2(pc(:,1:n_arc-1),arc);
dists=dists';

for i=1:length(analytes)
    for j=1:n_arc
       [tmp, p]=corrcoef(dists(:,j), DataPoints(:,i));
       corr_coef(i,j)=tmp(1,2);
       P_coef(i,j)=p(1,2);
    end
end

[sort_pc1, pc1_idx] = sort(corr_coef(:,1), 'descend');
[sort_pc2, pc2_idx] = sort(corr_coef(:,2), 'descend');
[sort_pc3, pc3_idx] = sort(corr_coef(:,3), 'descend');
[sort_pc4, pc4_idx] = sort(corr_coef(:,4), 'descend');

names_1=analytes(pc1_idx);
names_2=analytes(pc2_idx);
names_3=analytes(pc3_idx);
names_4=analytes(pc4_idx);

P_1=P_coef(:,1); P_1=P_1(pc1_idx);
P_2=P_coef(:,2); P_2=P_2(pc2_idx);
P_3=P_coef(:,3); P_3=P_3(pc3_idx);
P_4=P_coef(:,4); P_4=P_4(pc4_idx);


figure; 
barh(sort_pc1, 'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
hold on;
x=corr_coef(:,4);
barh(x(pc1_idx), 'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);
set(gca, 'YTick', 1:67, 'YTickLabel', names_1);


figure; 
barh(sort_pc2, 'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
x=corr_coef(:,3);
hold on;
barh(x(pc2_idx),'FaceColor',[0 0.447058826684952 0.74117648601532]);
set(gca, 'YTick', 1:67, 'YTickLabel', names_2);

figure; 
barh(sort_pc4,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);
hold on;
x=corr_coef(:,3);
barh(x(pc4_idx),'FaceColor',[0 0.447058826684952 0.74117648601532]);
set(gca, 'YTick', 1:67, 'YTickLabel', names_4);


% figure S7
coef_1=corr_coef(:,1); 
coef_2=corr_coef(:,2); 
coef_3=corr_coef(:,3);  
coef_4=corr_coef(:,4); 

sort_like_1=[sort_pc1, coef_2(pc1_idx), coef_3(pc1_idx), coef_4(pc1_idx),];
%figure; imagesc(sort_like_1);

%figure; imagesc(sort_like_1(:, [1,4]));
figure; imagesc(sort_like_1([1:10, 58:67], [1,4]));
set(gca, 'YTick', 1:20, 'YTickLabel',  fliplr(names_1([1:10, 58:67])));
colormap(flipud(spring))


sort_like2=[sort_pc2, coef_3(pc2_idx)];
figure; imagesc(flipud(sort_like2([1:10, 58:67],:)));
set(gca, 'YTick', 1:20, 'YTickLabel', fliplr(names_2([1:10, 58:67])));
colormap(flipud(spring))


sort_like4=[sort_pc4, coef_3(pc4_idx)];
figure; imagesc(flipud(sort_like4([1:10, 58:67],:)));
set(gca, 'YTick', 1:20, 'YTickLabel', fliplr(names_4([1:10, 58:67])));
colormap(flipud(spring))


figure; imagesc(flipud(sort_like_1([1:10, 58:67], :)));
set(gca, 'YTick', 1:20, 'YTickLabel', fliplr(names_1([1:10, 58:67])));
colormap(flipud(spring))



% figure S6
figure;  

[arc1_s, ia1] =sort(arcOrig(1,:), 'descend');
[arc2_s, ia2] =sort(arcOrig(2,:), 'descend');
[arc3_s, ia3] =sort(arcOrig(3,:), 'descend');
[arc4_s, ia4] =sort(arcOrig(4,:), 'descend');

subplot(1,4,1); barh(fliplr(arc1_s), ...
    'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
set(gca, 'YTick', 1:67, 'YTickLabel', analytes(ia1));
subplot(1,4,2); barh(fliplr(arc2_s),...
     'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(gca, 'YTick', 1:67, 'YTickLabel', analytes(ia2));
subplot(1,4,3); barh(fliplr(arc3_s), ...
    'FaceColor',[0 0.447058826684952 0.74117648601532]);
set(gca, 'YTick', 1:67, 'YTickLabel', analytes(ia3));
subplot(1,4,4); barh(fliplr(arc4_s), ...
    'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);
set(gca, 'YTick', 1:67, 'YTickLabel', analytes(ia4));

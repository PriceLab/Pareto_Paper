clear

load('IDs.mat');
load('Data.mat');
load('Parti_vars.mat');
load('rand_vec_7.mat');

data_not_norm=data(rand_vec, :);
mean_data=nanmean(data_not_norm);
std_data=nanstd(data_not_norm);

for i=1:size(data,2)
    col=data(:,i);
    col(isnan(col))=mean_data(i);
    data(:,i)=col;
end

[coefs1,scores1,variances] = pca(DataPoints);

ids=IDs;
visits=zeros(length(ids), 1);
i=2;
v=1;

while i<=length(ids)
    visits(i-1)=v;
    if ids(i)==ids(i-1)
        v=v+1;
    else
        v=1;
    end
    i=i+1;
end
if ~visits(end)
    visits(end)=1;
end

% figure S13 The distribution of trajectory length
visit_hist=[sum(visits==1), sum(visits==2), sum(visits==3), sum(visits==4), sum(visits==5), sum(visits==6),sum(visits==7), sum(visits==8),sum(visits==9)];
figure; 
bar(visit_hist); 
ylabel({'participants'});
xlabel({'Timepoints'});
title({'The distribution of trajectory length'});
set(gca,'FontSize',16,'XGrid','on','XTick',[1 2 3 4 5 6 7 8],'YGrid','on');

% dist between every two timepoint
dists=[];
for i=2:8
    indx=find(visits==i);
    for j=1:length(indx)
        data_i=data((indx(j)-1):indx(j), :);
        data_i=(data_i-repmat(mean_data, size(data_i,1),1))./repmat(std_data, size(data_i,1),1);
        data_i=data_i*coefs1(:,1:3);
        
        dists(end+1)=sqrt(sum((data_i(1,:)-data_i(2,:)).^2)); 
    end
end

% Fig S9: The distribution of distances between two timepoints
figure; 
hist(dists,30);
title({'Distribution of distances between visits'});
% [mean(dists), median(dists), std(dists), min(dists), max(dists)]



dists_fl=[];
long_dists=[];
long_dists_i=1;
for i_visits=3:8
    
    v_ind=find(visits==i_visits);
    
    for i=1:length(v_ind)
        indx=(v_ind(i)-(i_visits-1)):v_ind(i);
        data_i=data(indx, :);
        data_i=(data_i-repmat(mean_data, size(data_i,1),1))./repmat(std_data, size(data_i,1),1);
        data_i(isnan(data_i))=0;
        
        data_i=data_i*coefs1(:,1:3);
        
            tmp=sqrt(sum((data_i(1,:)-data_i(end,:)).^2));
            dists_fl(end+1)=tmp;
            long_dists(long_dists_i,1)=i_visits;
            long_dists(long_dists_i,2)=v_ind(i);
            long_dists(long_dists_i,3)=ids(v_ind(i));
            long_dists(long_dists_i,4)=tmp;
            long_dists_i=long_dists_i+1;
    end

end

figure; hist(long_dists(:,4),30);
title({'Distribution of distances v1-vn'});
%[mean(long_dists(:,4)), median(long_dists(:,4)), std(long_dists(:,4)), min(long_dists(:,4)), max(long_dists(:,4))]

%% stats on movements 

[max_visit, i_long_dists]=unique(long_dists(:,3),'last'); % choose the longest trajectory per ID
% there are 1186 unique trajectories of 3 or more timepoints 

long_dists_u=long_dists(i_long_dists,:);
[sorted_ld, sorted_ld_i]=sort(long_dists_u(:,4), 'descend'); % rank the trajecories according to delta firstV, lastV
long_dists_u=long_dists_u(sorted_ld_i,:);

% calculate the direction and consistancy of the movements

distFromArc1=ones(length(long_dists_u), 8).*NaN;
distFromArc2=ones(length(long_dists_u), 8).*NaN;
distFromArc3=ones(length(long_dists_u), 8).*NaN;
distFromArc4=ones(length(long_dists_u), 8).*NaN;


for i=1:length(long_dists_u)
    
    eve_idx=find(ids==long_dists_u(i,3));

    curr_data=data(eve_idx, :);
   
    curr_data(isnan(curr_data))=0;
    curr_data_norm=(curr_data-repmat(mean_data, size(curr_data,1),1))./repmat(std_data, size(curr_data,1),1);

    pos_on_tetra=curr_data_norm*coefs1(:,1:3);
    
    [~,tmp]=sortDataByDistance2(pos_on_tetra,arc);
    
    distFromArc1(i, 1:long_dists_u(i,1))=tmp(1,:);
    distFromArc2(i, 1:long_dists_u(i,1))=tmp(2,:);
    distFromArc3(i, 1:long_dists_u(i,1))=tmp(3,:);
    distFromArc4(i, 1:long_dists_u(i,1))=tmp(4,:);

end

delta1=[(distFromArc1(:,1)-distFromArc1(:,2)), ...
    (distFromArc1(:,2)-distFromArc1(:,3)), (distFromArc1(:,3)-distFromArc1(:,4)), ... 
    (distFromArc1(:,4)-distFromArc1(:,5)), (distFromArc1(:,5)-distFromArc1(:,6)), ...
    (distFromArc1(:,6)-distFromArc1(:,7)), (distFromArc1(:,7)-distFromArc1(:,8))];

delta2=[(distFromArc2(:,1)-distFromArc2(:,2)), ...
    (distFromArc2(:,2)-distFromArc2(:,3)), (distFromArc2(:,3)-distFromArc2(:,4)), ... 
    (distFromArc2(:,4)-distFromArc2(:,5)), (distFromArc2(:,5)-distFromArc2(:,6)), ...
    (distFromArc2(:,6)-distFromArc2(:,7)), (distFromArc2(:,7)-distFromArc2(:,8))];

delta3=[(distFromArc3(:,1)-distFromArc3(:,2)), ...
    (distFromArc3(:,2)-distFromArc3(:,3)), (distFromArc3(:,3)-distFromArc3(:,4)), ... 
    (distFromArc3(:,4)-distFromArc3(:,5)), (distFromArc3(:,5)-distFromArc3(:,6)), ...
    (distFromArc3(:,6)-distFromArc3(:,7)), (distFromArc3(:,7)-distFromArc3(:,8))];

delta4=[(distFromArc4(:,1)-distFromArc4(:,2)), ...
    (distFromArc4(:,2)-distFromArc4(:,3)), (distFromArc4(:,3)-distFromArc4(:,4)), ... 
    (distFromArc4(:,4)-distFromArc4(:,5)), (distFromArc4(:,5)-distFromArc4(:,6)), ...
    (distFromArc4(:,6)-distFromArc4(:,7)), (distFromArc4(:,7)-distFromArc4(:,8))];

delta1_isgrater0=sum(delta1>0,2);
delta2_isgrater0=sum(delta2>0,2);
delta3_isgrater0=sum(delta3>0,2);
delta4_isgrater0=sum(delta4>0,2);

% this analysis shows that participants move to any direction, but more
% away from arc4 - the unhealthy archetype, and only for small number of 
% timepoint away from archetype 1 - the older and the healthy
[sum(delta1_isgrater0==0), sum(delta2_isgrater0==0), sum(delta3_isgrater0==0), ...
    sum(delta4_isgrater0==0)]

[sum(((long_dists_u(:,1)-1)-delta1_isgrater0)==0), sum(((long_dists_u(:,1)-1)-delta2_isgrater0)==0), ...
    sum(((long_dists_u(:,1)-1)-delta3_isgrater0)==0), sum(((long_dists_u(:,1)-1)-delta4_isgrater0)==0)]



% out of the top 2% of longest trajectories only 3 are moving away from
% archetype 1, and 21 are moving away from archetype 4 as expected from a
% wellness program. The 3 unique trajectories are described in the MS
away_from_arc1=zeros(24,1);
away_from_arc4=zeros(24,1);

for i=1:24
    n_v=long_dists_u(i,1);
    if (distFromArc4(i,1)-distFromArc4(i,n_v))<0
       away_from_arc4(i)=1;
    end
    
    if (distFromArc1(i,1)-distFromArc1(i,n_v))<0
       away_from_arc1(i)=1;
    end
end

% [sum(away_from_arc1), sum(away_from_arc4)]


% 25/1186 trajectories had a time point that exceeded the tetrahedron
% boundaries. There was one trajectory that had 3 timpeoints, and one that
% had 5/6 timpoints out of the convex hull. These trajectories are
% discussed in the paper

% count how many time points exceeded the polygon 
count_exceed=zeros(length(long_dists_u),1);

for i=1:length(long_dists_u)
    idx=find(ids==long_dists_u(i,3));

    curr_data=data(idx, :);
   
    curr_data(isnan(curr_data))=0;
    curr_data_norm=(curr_data-repmat(mean_data, size(curr_data,1),1))./repmat(std_data, size(curr_data,1),1);

    pos_on_tetra=curr_data_norm*coefs1(:,1:3);
    count_exceed(i)= sum(~inhull(pos_on_tetra, arc));
end


%[sum(count_exceed>0), sum(count_exceed==3), sum(count_exceed==4), sum(count_exceed==5)]



% Fig 6: plotting the unique trajectories on the tetrahedron

close all
pos=[1:3, 6, find(count_exceed==5)];

for i=1:5
    
    rows=long_dists_u(pos(i),2)-long_dists_u(pos(i),1)+1:long_dists_u(pos(i),2);
    curr_data=data(rows,:);
    
    curr_data_norm=(curr_data-repmat(mean_data, size(curr_data,1),1))./repmat(std_data, size(curr_data,1),1);
    
    pos_on_tetra=curr_data_norm*coefs1(:,1:3);
    
    uiopen('empty_tetra.fig',1)
    hold on
    plot3(pos_on_tetra(:,1),pos_on_tetra(:,2),pos_on_tetra(:,3),'*-m');
    hold on;
    plot3(pos_on_tetra(1,1),pos_on_tetra(1,2),pos_on_tetra(1,3),'ok');
    
    
end


% generate Figure S11
close all

load('clinical_labs_ls.mat')
i=find(count_exceed==5);
rows=long_dists_u(i,2)-long_dists_u(i,1)+1:long_dists_u(i,2);
curr_data=data(rows,:);
curr_data_norm=(curr_data-repmat(mean_data, size(curr_data,1),1))./repmat(std_data, size(curr_data,1),1);

figure; 
subplot(3,1,1); boxplot(DataPoints(:, 1:23),analytes(1:23), 'Symbol', '.k','OutlierSize',2);
xtickangle(45);
ylim([-10 15]);
hold on; plot(curr_data_norm(2:end ,1:23)', 'oc');

subplot(3,1,2);
boxplot(DataPoints(:, 24:46),analytes(24:46), 'Symbol', '.k','OutlierSize',2);
xtickangle(45);
ylim([-10 15]);
hold on; plot(curr_data_norm(2:end ,24:46)', 'oc');

subplot(3,1,3);
boxplot(DataPoints(:, 47:67),analytes(47:67), 'Symbol', '.k','OutlierSize',2);
xtickangle(45);
ylim([-10 15]);
hold on; plot(curr_data_norm(2:end ,47:67)', 'oc');

% generate Figure S12
close all

load('clinical_labs_ls.mat')
i=2;
rows=long_dists_u(i,2)-long_dists_u(i,1)+1:long_dists_u(i,2);
curr_data=data(rows,:);
curr_data_norm=(curr_data-repmat(mean_data, size(curr_data,1),1))./repmat(std_data, size(curr_data,1),1);

figure; 
subplot(3,1,1); boxplot(DataPoints(:, 1:23),analytes(1:23), 'Symbol', '.k','OutlierSize',2);
xtickangle(45);
ylim([-10 15]);
hold on; plot(curr_data_norm(2:end ,1:23)', 'oc');

subplot(3,1,2);
boxplot(DataPoints(:, 24:46),analytes(24:46), 'Symbol', '.k','OutlierSize',2);
xtickangle(45);
ylim([-10 15]);
hold on; plot(curr_data_norm(2:end ,24:46)', 'oc');

subplot(3,1,3);
boxplot(DataPoints(:, 47:67),analytes(47:67), 'Symbol', '.k','OutlierSize',2);
xtickangle(45);
ylim([-10 15]);
hold on; plot(curr_data_norm(2:end ,47:67)', 'oc');

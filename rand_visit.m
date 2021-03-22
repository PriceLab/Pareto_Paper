clear
load('IDs.mat');
load('Data.mat');

j=1;
new_db=[];
rand_vec=[];
i=1;
while i<length(IDs)
    idx=[];
    while i<length(IDs) && IDs(i)== IDs(i+1)
        idx=[idx, i];
        i=i+1;
    end
    
    if isempty(idx)
        new_db(j,:)=data_white(i,:);
        rand_vec(j)=i;
        j=j+1;
        i=i+1;
    else
        idx=[idx,i];
        pos = randi(length(idx));
        new_db(j,:)=data_white(idx(pos),:);
        rand_vec(j)=idx(pos);
        i=i+1;
        j=j+1;
    end
end


rand_vec(j)=7367;
new_db(j,:)=data_white(7367 ,:);

data_znorm=(new_db-repmat(nanmean(new_db), size(new_db,1), 1))./(repmat(nanstd(new_db), size(new_db,1), 1));

data_znorm_no_nans=data_znorm;
for i=1:size(data_znorm,2)
    col=data_znorm(:,i);
    col(isnan(col))=nanmean(col);
    data_znorm_no_nans(:,i)=col;
end

% [coeff,score,latent] = pca(data_znorm_no_nans);
% 
% figure; plot3(score(:, 1), score(:, 2),score(:, 3), 'o');
%

n_arcs=4;

[arc, arcOrig, pc] = ParTI_findArcs_tRatio(data_znorm_no_nans,1,10, 'w_alg1_fast_4A_randVx', n_arcs);

% save('ParTI_vars', 'data_znorm_no_nans', 'pc', 'arc', 'arcOrig');
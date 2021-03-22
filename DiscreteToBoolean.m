function [boolMat, featName] = DiscreteToBoolean(enriched, labels, cols, opf_ID)
%Enriched is a cell array. Each row is a sample, each column is a category
if ~iscell(enriched)
    enriched = cellfun(@(x) num2str(x), (num2cell(enriched)),'UniformOutput', false);
end

if length(labels)~=size(enriched,2)
   boolMat = NaN;
   featName = NaN;
   fprintf(opf_ID, 'Error - number of columns of Enriched should match number of category-names');
   return;
end
  
if nargin == 2 || cols==0
    cols = 1:size(enriched,2);
end
 enr = enriched(:,cols);
 lab = labels(cols);
 featName = {};
 boolMat = [];
for index = 1:size(enr,2)
    u = unique(enr(:,index));
     
    for subind = 1:length(u)    
        temp = ([lab{index},': ',u{subind}]);
        featName = [featName,temp];  
        boolMat =[boolMat,strcmp(enr(:,index),u{subind})];
    
    end
end

    
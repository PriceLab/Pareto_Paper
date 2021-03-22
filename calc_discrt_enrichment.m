function [] = calc_discrt_enrichment(data,Lia, features, fname)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
bin=data(Lia==1, :);
rest=data(Lia==0, :);


if sum(isnan(bin))==length(bin) 
    disp('error: no elements in bin');
    return
end



pval = 1 - hygecdf(sum(bin), length(Lia) , sum(data), sum(Lia));



table =  [pval', medianDifference', meanDifference', isSignificantAfterFDR', ...
          meanFirst', meanRest', medianFirst', medianRest'];
      
continuousTitles = { 'Feature Name', 'P value (Mann-Whitney)'...
,'Median Difference' ,'Mean Difference','Significant after Benjamini-Hochberg correction?',...
'Mean First',  'Mean Rest', 'Median First',  'Median Rest',}; 
[ordContTable, idx] = sortrows(table, 1 );
ordContFeatureNames = features(idx)';

fullTable=cell(size(rest,2)+1, length(continuousTitles));
fullTable(1,:)= continuousTitles; 
fullTable(2:end,1) = ordContFeatureNames;
fullTable(2:end, 2:end)=num2cell(ordContTable);

if nargin==4
    cell2csv(fname, fullTable);
end

end

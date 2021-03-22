function [ordering, distances]=sortDataByDistance(DataPoints,Archetypes)
%function [ordering, distances]=sortDataByDistance(DataPoints,Archetypes)
% 1. Data points is the values of different traits (e.g. expression
% level of genes) - each sample is a row, each trait (gene) is a column 
% 2. Archetypes is a matrix with different archetypes in rows

%initializing 
distances = zeros(size(Archetypes,1),size(DataPoints,1));
ordering = zeros(size(Archetypes,1),size(DataPoints,1));
dist= zeros(size(Archetypes,1),size(DataPoints,1));

for arch = 1:size(Archetypes)
    
    archAtOrigin = bsxfun(@minus,DataPoints,Archetypes(arch,:));
    % archAtOrigin moves the archetype to the origin
    dist(arch,:) = sqrt(sum(archAtOrigin.^2,2));
    % dist is the Euclidian distance between the datapoints to
    % archetype(arch);
    [sortedDistances,pos] = sort( dist(arch,:));
    % the distances with the corresponding indexes
    distances(arch,:) = sortedDistances';
    ordering(arch,:) = pos';
    
%     distances=dist;
%     ordering=[];
end

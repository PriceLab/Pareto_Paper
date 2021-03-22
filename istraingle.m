function pval = istraingle(smallmatpc,T,pruns,delta)
% T - t ratio to compare with
if nargin<4
    delta = 0;
end

T_shuf = zeros(pruns,1);
arcnum = size(smallmatpc,2)+1;
for j=1:pruns
    % shuffle data
    newmat = ShuffleData(smallmatpc);
    % compute archetypes
    U=1:size(newmat,1); % Entries in X used that is modelled by the AA model
    I=1:size(newmat,1); % Entries in X used to define archetypes
    [arc,S,C,SSE,varexpl]=PCHA1(newmat',arcnum,I,U,delta);
    % compute T ratio
     T_shuf(j) = t_ratio(newmat,arc); 
     if ~mod(j,10)
     j
     end

end

pval = sum(T_shuf>=T)/length(T_shuf);


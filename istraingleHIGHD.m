function [pval,effectsize] = istraingleHIGHD(data,arcnum,pruns,delta)
% data - samples X coordinates
% arcnum - number of archetypes
% pruns - number of shuffled data sets
% delta - PCHA parameter

if nargin<4
    delta = 0;
end

% Compute real EV
U=1:size(data,1); % Entries in X used that is modelled by the AA model
I=1:size(data,1); % Entries in X used to define archetypes
[~,~,~,~,varexplR]=PCHA1(data',arcnum,I,U,delta);

% Compute random EV vector
EV_shuf = zeros(pruns,1);
for j=1:pruns
    % shuffle data
    newmat = ShuffleData(data);
    % compute shuffled EV
    U=1:size(newmat,1); % Entries in X used that is modelled by the AA model
    I=1:size(newmat,1); % Entries in X used to define archetypes
    [~,~,~,~,varexpl(j)]=PCHA1(newmat',arcnum,I,U,delta);
     if ~mod(j,10)
     j
     end
end

pval = sum(varexpl >= varexplR)/length(varexpl);
effectsize = (varexplR-mean(varexpl))/std(varexpl) ;


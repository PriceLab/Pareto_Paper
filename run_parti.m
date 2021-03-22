
clear
ls=dir('Parti_data/*.mat');
algNum=1;
dim=10;
cols=0;
binSize=0.05;
n_arc=4;

 % To save computational time, after the archetype position was determined
 % their position and and all relevant vars were saved and uploaded for
 % further analysis. to calculate the archetype position run ParTI.m with
 % the proper variables, and uncomment line 123
 

f_vars='Parti_vars.mat';
load(f_vars);


for i=1:length(ls)
    load(['parti_data/',ls(i).name]);
    disp(ls(i).name);
    fname=regexprep(ls(i).name, '.mat', '');
    OutputFileName=['Parti_results/', fname ];
    
    ParTI(DataPoints,algNum,dim,DiscFeatName,EnMatDis,cols,ContFeatName,EnMatCont,[],binSize,OutputFileName, n_arc, f_vars);
end






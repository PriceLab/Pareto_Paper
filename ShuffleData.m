function newmat = ShuffleData(data)
    newmat = zeros(size(data));  
    for i=1:size(data,2)
        newmat(:,i) = randsample(data(:,i),size(data,1),true); %sampling with replacement
    end
end
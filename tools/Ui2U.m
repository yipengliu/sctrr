function U = Ui2U(Ui)
% Ui2U  Merge all the core factors in tensor ring regression into its original
% tensor.
%  Input:
%       Ui: the core factors in the TR decomposition
%  Output:
%       U: the origin tensor
%Tensor Ring Ridge Regression
%Copyright 2019
%

    d = length(Ui); % dimension    
    
    I = nan(1, d);      
    for i =1:d
        I(i) = size(Ui{i},2);% the size of each dimension
    end
    
    T = TCP(Ui(1:end-1));  % tensor connect product of the first d-1 core factors
    [r1,n,r2]=size(Ui{end});
    
    U= double(tenmat(T,2))*reshape(permute(Ui{end},[3,1,2]),r1*r2,n); % the relation between the TR decomposition and its origin tensor in Lemma 2.
    U = reshape(U, I); % reshape U into the origin size.
end
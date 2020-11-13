function X = project(U,lambda)
% function X = proj_l1(U, lambda)
% Description:
%   Soft Thresoding function. Solve one of the following problems:
%       U is a matrix of size d x k 
%       lambda is a positive scalar
%            X = argmin_X 0.5*||X - U|| + lambda||X||_1 
% Inputs: U: double dense matrix d x k 
%         lambda
% Outputs: X: a full matrix in d x k
% -----------------------------------------------

    %%
 if lambda==0
     X=U;
 else
     X = max(0, U - lambda) + min(0, U + lambda);
 end
 
end
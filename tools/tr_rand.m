function [tr]=tr_rand(n,d,r)
%Generates a random tensor
%   [Tr]=Tr_RAND(N,D,R) Generates a random tensor with dimensions specified
%   by (N,D), where N can be a number of an array of dimension D, R is a
%   rank or a array of dimension d+1
%---------------------------
if ( numel(n) == 1 ) 
 n = n * ones(d,1);
end
if ( numel(r) == 1 )
 r = r * ones(d+1,1); 
end
U=cell(d,1);
U{1}=(randn(r(d+1),n(1),r(1)));
for i=2:d
U{i}=(randn(r(i-1),n(i),r(i)));
end
tr.n=n;
tr.r=r;
tr.d=d;
tr.U=U;
return
end

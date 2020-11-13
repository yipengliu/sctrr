function T=self_contract(U)
% self_contract is an extension of U'U in matrix space into tensor space. 
% It contracts U and itself along the middle index.
%  input:  U in the size of [l,c,r]
%  output: T in the size of [l*l,r*r]

[l, c, r] = size(U);
T = reshape(permute(U,[2,1,3]), [c, l*r]);
T=T'*T;
T=permute(reshape(T,[l,r,l,r]),[1,3,2,4]);
T=reshape(T,[l*l,r*r]);
end
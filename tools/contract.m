function[Y]=contract(X,W,L)
% compute the contract product between two tensors.
%    Y= X * W

dimx=size(X);
N=dimx(1);
P=dimx(2:L+1);
dimmodel=size(W);
Q=dimmodel(L+1:length(dimmodel));
Xmat=reshape(X,[N,prod(P)]);
Wmat=reshape(W,[prod(P),prod(Q)]);
Ymat=Xmat*Wmat;
dimy=[N,Q];
Y=reshape(Ymat,dimy);

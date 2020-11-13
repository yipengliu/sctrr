function Ul=update_ul(X,Us,Ub,V,Ymat,Lam,lambda,index)

 s=size(Ymat,1);
[rl,q,rr]=size(V);
Us=permute(Us,[1,3,2]);
Ub=permute(Ub,[2,1,3]);
% COMPUTE computational complexity of diffrent contraction order, choose
% the smallest one.
t1=(rl*rr)^2*(s*prod(index)+q)+(rl*rr)*s*prod(index)^2;
t2=(rl*rr)^2*((s+1)*prod(index)^2+q);
t3=(s*prod(index)*q *(rl*rr*+prod(index)));
t=min(min(t1,t2),t3);
switch t
    case  t1    
        Vs=double(tenmat(V,2));
        Yred=reshape((Ymat*Vs),[s,rl,rr]);
        A=TCP({Us,X,Ub});
        A=reshape(A,[rr,index(2),s,index(1),index(3),rl]);
        A=permute(A,[4,2,5,3,6,1]);
        A=reshape(A,[prod(index),s,rl*rr]);
        temp1=reshape(A,prod(index),[])*Yred(:); % COMPUTE ATY
        ATA=reshape(A,[],rl*rr)*(Vs'*Vs);
        ATA=reshape(ATA,[prod(index),s,rl*rr]);
        ATA=reshape(A,[],s*rl*rr)*reshape(ATA,prod(index),[])';
        Ul=(ATA+lambda.*Lam)\(temp1);
        Ul=permute(reshape(Ul,index),[2,1,3]);
    case t2      
        Vs=double(tenmat(V,2));
        Yred=reshape((Ymat*Vs),[s,rl,rr]);
        A=TCP({Us,X,Ub});
        A=reshape(A,[rr,index(2),s,index(1),index(3),rl]);
        A=permute(A,[4,2,5,3,6,1]);
        A=reshape(A,[prod(index),s,rl*rr]);
        temp1=reshape(A,prod(index),[])*Yred(:);% COMPUTE ATY
        ATA=compute_UiTUi({A,Vs'});
        Ul=(ATA+lambda.*Lam)\(temp1);
        Ul=permute(reshape(Ul,index),[2,1,3]);
    case t3   
        A=Ui2U({Us,X,Ub,V});
        A=reshape(A,[index(2),s,index(1),index(3),q]);
        A=permute(A,[2,5,3,1,4]);
        A=reshape(A,s,q,[]);
        ATA=reshape(A,s*q,[]);
        ATA=ATA'*ATA;
        Ul=(ATA+lambda.*Lam)\(A'*Ymat(:));
        Ul=permute(reshape(Ul,index),[2,1,3]);
end
end
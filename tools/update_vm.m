function Vm=update_vm(Y,Bcell,Lam,lambda)

   [r1,n,r2]=size(Bcell{1});% r1 * (I2...In) * rn
   
    B = TCP(Bcell(2:end)); 
    B=reshape(permute(B,[2,3,1]),[],r1*r2);
    BTB=compute_UiTUi(Bcell(2:end));% fast contraction of B'B  
    
    
    Vm=(Y*B)/(BTB+lambda.*Lam);
    Vm=permute(reshape(Vm,[n,r1,r2]),[2,1,3]);
end
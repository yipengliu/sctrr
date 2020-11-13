%% Demo to test the proposed Smooth Compact Tensor Ring Regression
clear all;

d=5;L=3;M=2;
N=1000;
nlist=[5,10,15,20];
Rlist={[2,4,5,3,4,4],repmat(3,[1,d+1]),repmat(5,[1,d+1])};
for i=1
    n=nlist(i);
    P=repmat(n,[1,L]);
    Q=repmat(n,[1,M]);
    dim=[P,Q];
    for j=1
        r=Rlist{j};
        %% the model
        model = tr_rand(dim,d,r);
        om=Ui2U(model.U);
        
        X=randn([N,P]);
        Y=contract(X,om,3)+1*random('Normal', 0, 1, [N,Q]);
        
        XS=randn([N,P]);
        YS=contract(XS,om,3)+1*random('Normal', 0, 1, [N,Q]);
        
        para.maxiter =100;
        
        %         tic;
        
        para.lambda=10000;
        MaxRank=10;
        para.tol=1e-3;
        beta=0.1;
        [model,runtime] = SCTRR(para, X,Y,MaxRank,beta);
        w=Ui2U(model.U);
        
        %         toc;
        estimated_model_error(i,j)=norm(om(:)-w(:),'fro')/norm(om(:),'fro');
        Ypred=contract(XS,w,3);
        Ypred=(reshape(Ypred,[numel(YS),1]));
        Y=(reshape(YS,[numel(YS),1]));
        cor(i,j)= mycorrcoef(Ypred(:),Y(:));
        Ypress = sum((Y(:)-Ypred(:)).^2);
        rmse(i,j) = sqrt(Ypress./numel(Y));
        Q2(i,j) = 1 - Ypress./sum(Y(:).^2);
        
    end
end
function [model, runtime]=SCTRR(para, X,Y,MaxRank,scale)
%

%% set parameters
MaxAlpha = 1e8;
maxiter = para.maxiter;
Update_Method='grl';
lambda=para.lambda;


sx=size(X);
sy=size(Y);
if sx(1)==sy(1)
    N=sx(1);
else
    print('error: the sample number is not the same');
    return;
end
P=sx(2:end);
Q=sy(2:end);
L=length(P);
M=length(Q);
dim=[P,Q];
d=L+M;
MinU =para.tol;

AlphaEps = 1e-3;
noiterPruning = 1; % the iteration to begin with prunning the component;


%% initialize parameters
model = init_trmodel(dim,MaxRank.*ones(d+1,1)); % the latent factors for the coefficent tensor


Alpha=cell(d,1);Keep_indx=cell(d,1);
for i=1:d
    Alpha{i}= ones(MaxRank,1);
    Keep_indx{i} = [1:MaxRank]';
end


tic
Xmat = reshape(X,[N,prod(P)]);
Ymat = reshape(Y,[N,prod(Q)]);



%% compute Xarray
Xarray=cell(1,L);
for i=1:L
    ps=prod(P(1:i-1));
    pb=prod(P(i+1:L));
    Xarray{i}=reshape(permute(reshape(X,[N,ps,P(i),pb]),[2,1,3,4]),[ps,N*P(i),pb]);
end

%% compute Us, Ub
r=model.r;
Us=cell(1,L);
Ub=cell(1,L);
Us{1}=reshape(eye(r(d+1),r(d+1)),[r(d+1),1,r(d+1)]);
Ub{L}=reshape(eye(r(L),r(L)),[r(L),1,r(L)]);
for i=2:L
    Us{i}=reshape(left_unfold(Us{i-1})*right_unfold(model.U{i-1}),[r(d+1),prod(P(1:i-1)),r(i-1)]);
end
for i=L-1:-1:1
    Ub{i}=reshape(left_unfold(model.U{i+1})*right_unfold(Ub{i+1}),[r(i),prod(P(i+1:L)),r(L)]);
end

%% Main LOOP

for k=1:maxiter
    %% adjust the rank 
    if k > noiterPruning
        r_old=r;
        for i=1:d
            if r(i) ~= 1
                if max(Alpha{i}) > MaxAlpha-AlphaEps
                    L_index = Alpha{i} <= MaxAlpha-AlphaEps;
                    Keep_indx{i} = Keep_indx{i}(L_index);
                    r(i) = sum(L_index==1);
                    if r(i) == 0
                        fprintf('No efficient rank, Pleae Inrease the Max Rank');
                        break;
                    end
                    Um = model.U{i};
                    Um(:,:,~L_index) = [];
                    model.U{i}= Um;
                    if i==d
                        Um = model.U{1};
                        Um(~L_index,:,:) = [];
                        model.U{1}= Um;
                    else
                        Um = model.U{i+1};
                        Um(~L_index,:,:) = [];
                        model.U{i+1}= Um;
                    end
                    Alpha{i} = Alpha{i}(L_index);
                end
            end
        end
        r(d+1)=r(d);
        if any(r<r_old)
            Us{1}=reshape(eye(r(d+1),r(d+1)),[r(d+1),1,r(d+1)]);
            Ub{L}=reshape(eye(r(L),r(L)),[r(L),1,r(L)]);
            for i=2:L
                Us{i}=reshape(left_unfold(Us{i-1})*right_unfold(model.U{i-1}),[r(d+1),prod(P(1:i-1)),r(i-1)]);
            end
            for i=L-1:-1:1
                Ub{i}=reshape(left_unfold(model.U{i+1})*right_unfold(Ub{i+1}),[r(i),prod(P(i+1:L)),r(L)]);
            end
            
        end
    end
    
    %% update U1...Ud
    U_old=reshape(left_unfold(Us{L})*right_unfold(model.U{L}),[r(d+1),prod(P(1:L)),r(L)]);
    Vs=TCP(model.U(L+1:L+M));
    for i=1:L
        [r1,n,r2]=size(model.U{i});
        if i==1
            Lam{1} = diag(Alpha{d}./scale');
            Lam{2} = diag(Alpha{1}./scale');
        else
            Lam{1} = diag(Alpha{i-1}./scale');
            Lam{2} = diag(Alpha{i}./scale');
        end
        Lam_total=kron((kron(Lam{2}',eye(size(Lam{1})))+kron(eye(size(Lam{2})),Lam{1})),eye(n,n));
        model.U{i}=update_ul(Xarray{i},Us{i},Ub{i},Vs,Ymat,Lam_total,lambda,[n,r1,r2]);
        if i~=L
            Us{i+1}=reshape(left_unfold(Us{i})*right_unfold(model.U{i}),[r(d+1),prod(P(1:i)),r(i)]);
        end
        if i~=1
            Ub{i-1}=reshape(left_unfold(model.U{i})*right_unfold(Ub{i}),[r(i-1),prod(P(i:L)),r(L)]);
        end
    end
    
    U = TCP(model.U(1:L));
    Xred=double(ttm(tensor(U),Xmat,2));
    Bcell{1,1}=Xred;
    for i=1:M
        Lam{1} = diag(Alpha{i+L-1}./scale');
        Lam{2} = diag(Alpha{i+L}./scale');
        Bcell(2:M+1,1)=model.U(L+1:L+M);
        Lam_total=(kron(Lam{2}',eye(size(Lam{1})))+kron(eye(size(Lam{2})),Lam{1}));
        model.U{i+L}=update_vm(modelk_unfold(Y,i+1),Bcell([i+1:M+1,1:i]),Lam_total,lambda);
    end
    
    
    %% update parameter Alpha
    if k > noiterPruning
        for i=1:d
            if i==d
                UA = sum(left_unfold(model.U{d}).^2,1) + sum(right_unfold(model.U{1})'.^2,1);
            else
                UA = sum(left_unfold(model.U{i}).^2,1) + sum(right_unfold(model.U{i+1})'.^2,1);
            end
            
            
            if strcmp(Update_Method, 'grl')
                alpha = sqrt(UA);
                alpha = alpha./sum(alpha) + 1e-12;
                Alpha{i} = 1./alpha;
            elseif strcmp(Update_Method, 'fix')
                Alpha{i} = ones(1, r);
            end
        end
    end
    %% Stop conditions
    % compute error
    
    V_new=TCP(model.U(L+1:L+M));
    
    NmU=(norm(U(:)-U_old(:))+norm(V_new(:)-Vs(:)))/(norm(U_old(:))+norm(Vs(:)));
    if isnan(NmU)
        break;
    end
    g=sprintf('%d ', r);
    fprintf('Iteration: %d, Weight Change: %0.5f, Rank: %s \n,', k , NmU,  g )
    if k>1 && ( NmU<MinU)
        break;
    end
    
end
model.r=r;
runtime=toc;
fprintf('running time=%fs\n',runtime);
end



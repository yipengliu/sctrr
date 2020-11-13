function T = compute_UiTUi(U)
    % compute_UiTUi caculates the network of
    %      
    %    r_0  ----U1----U2----...----Un---- r_n
    %             |     |      |     |
    %    r_0  ----U1----U2----...----Un---- r_n
    %
    %  first contract Ul and itself, then yields
    %
    %    r_0*r_0  ----U1'U1----U2'U2----...----Un'Un---- r_n*r_n
    %
    % Input: U is a cell {U1, U2,..., Un}
    %        U_i is a tensor:       r_{i-1} * I_i * r_i
    % Output: T is a tensor:         [r_0 * r_n,r_0 * r_n]
    
    n = length(U);
    [l,~,~]=size(U{1});
    [~,~,r]=size(U{n});
    for i=1:n
    U{i} = self_contract(U{i});
    end
      T=U{1};
    if(n==1)  
         T=reshape(permute(reshape(T,[l,l,r,r]),[3,1,4,2]),[l*r,l*r]);
        return;
    end
    
    for i =2 : n
        T = T*U{i};
    end
    T=reshape(permute(reshape(T,[l,l,r,r]),[3,1,4,2]),[l*r,l*r]);
end
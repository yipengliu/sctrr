function MU = middle_unfold(U,flag)
    % input:    U is a 3rd order tensor     r_i-1 * I_i * r_i
    % output:   LU is a matrix              r_i-1 I_i * r_i
    if flag=='p'
    MU = reshape(permute(U,[2,3,1]), size(U,2),[]);
    elseif flag=='r'
    MU = reshape(permute(U,[2,1,3]), size(U,2),[]);   
    else
     return      
    end
end
function mat= modelk_unfold(U,k)
    % input:    U is a 3rd order tensor     r_i-1 * I_i * r_i
    % output:   RU is a matrix              r_i-1 * I_i r_i
    
    d=length(U);
    n=size(U,k);

    mat = reshape(permute(U,[k:d,1:k-1]),n,[]);
end
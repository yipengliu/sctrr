function mat= L_unfold(U,L)
    % input:    U is a 3rd order tensor     r_i-1 * I_i * r_i
    % output:   RU is a matrix              r_i-1 * I_i r_i
    dim=size(U);
    mat = reshape(U, prod(dim(1:L)),[]);
end
function mat= model1_unfold(U)
    % input:    U is a 3rd order tensor     r_i-1 * I_i * r_i
    % output:   RU is a matrix              r_i-1 * I_i r_i
  
    mat = reshape(U, size(U,1),[]);
end
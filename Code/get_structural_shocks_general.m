function [structural_shocks, IR] = get_structural_shocks_general(A,gamma,resid,which_shocks)

D = [gamma null(gamma')]; % building D, the only orthogonal matrix
IR = A*D; 
s = IR^(-1)*resid';
s = s'; % to make it (T,nshocks_recovered)
structural_shocks = s(:,1:length(which_shocks));
end
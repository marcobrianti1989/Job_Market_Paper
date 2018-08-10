  function [c,ceq] = constraint_orthogonality(gam)
       
         c = [];
         ceq = gam'*gam-eye(size(gam,2));
         
   end
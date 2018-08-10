   function [c,ceq] = constraint_barskysims(gam) 
         c = [];
         ceq = gam'*gam - 1;   
         
   end
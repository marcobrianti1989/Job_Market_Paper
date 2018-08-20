sig1 = sigma1(1,1);
sig12 = sigma1(2,1);
sig2 = sigma1(2,2);


A11 = sig1^0.5;
A21 = sig12/(sig1^0.5);
A22 = (sig2 - sig12^2/sig1)^0.5;

A = [A11 0; A21 A22]




























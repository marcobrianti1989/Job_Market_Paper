function [obj] = calibratef(params,xt,stdterm,meanterm)
% -------------------------------------------------------------------------
% Calibrate principal component to have mean = meanterm, std = stdterm
% -------------------------------------------------------------------------

a = params(1);
b = params(2);
term  =  exp((a.*xt+b)./2);
obj1  = (std(term)-stdterm)^2;
obj2  = (mean(term)-meanterm)^2;
obj   = obj1+obj2;
end
function [R2,mR2] = mrsq(Fhat,lamhat,ve2,colheaders,vartype,outfile)
% -------------------------------------------------------------------------
% Compute R2 and marginal R2 from estimated factors and factor loadings
%   Input:  Fhat       = estimated factor matrix
%           lamhat     = factor loadings
%           ve2        = eigenvalues
%           colheaders = variable names 
%           vartype    = variable type
%           outfile    = name of output file
%   Output: R2         = R2 for each series
%           mR2        = marginal R2 for each group
% -------------------------------------------------------------------------

% Compute values
[N ic]     = size(lamhat);                          
R2         = zeros(N,ic);                           
mR2        = zeros(N,ic);
tcode      = vartype;
for i = 1:ic
    R2(:,i)  = (var(Fhat(:,1:i)*lamhat(:,1:i)'))';  
    mR2(:,i) = (var(Fhat(:,i:i)*lamhat(:,i:i)'))';
end

% Print to file
fp          = fopen(outfile,'w');
in.fid      = fp;
in.outwidth = N;
group_id    = [ones(2,1);ones(3,1)*4;ones(15,1);...    %group identifiers
               ones(30,1)*2;ones(10,1)*3;...
               ones(10,1)*4;ones(11,1)*5; ...
               ones(4,1)*8; ones(22,1)*6;...
               ones(21,1)*7;ones(3,1)*2;4;...
               ones(N-132,1)*9];
series      = [];
for i = 1:N
    series  = str2mat(series,colheaders{i});
end
series      = series(2:end,:);
fmt1        = ['\n                            ' repmat('%7d',1,ic)];
fmt         = ['\n' repmat('%4d',1,3)  ' %s' repmat('%7.3f',1,ic) ];
fprintf(fp,fmt1,(1:1:ic));

fprintf(fp,'\n R2\n');
for i = 1:N;
    fprintf(fp,fmt, i,group_id(i),tcode(i),series(i,:),R2(i,1:ic));
end;
fprintf(fp,'\n\n Marginal R2\n');
for i = 1:N;
    fprintf(fp,fmt,i,group_id(i),tcode(i),series(i,:),mR2(i,1:ic));
end;

[vals,ind] = sort(mR2,'descend');
Rsq        = ve2./sum(ve2);
fprintf(fp,'\n\n Highest Marginal R2''s');
for ii = 1:ic;
    fprintf(fp,'\n\n Factor %d: total contribution = %0.4f \n',ii,Rsq(ii));
    for j = 1:10;
        i = ind(j,ii);
        fprintf(fp,fmt,i,group_id(i),tcode(i),series(i,:),mR2(i,ii));
    end
end;
fclose(fp);
end
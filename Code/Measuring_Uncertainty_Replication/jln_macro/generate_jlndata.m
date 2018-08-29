% -------------------------------------------------------------------------
% Load raw data and generate cleaned data used for analysis
% -------------------------------------------------------------------------

% Load data
clear; clc;
% for octave
% pkg load io;
mc=1;
[mdata,mtxt] = xlsread('jlnrawdata.xlsx',1);
[fdata,ftxt] = xlsread('jlnrawdata.xlsx',2);
return
%%% windows users
windows=0; mac=1;
if windows==1;
vartype      = mdata(1,:);
mdata        = mdata(3:end,:);
mnames       = mtxt(2,2:end);
end;
if mac==1;
vartype      = mdata(1,2:end);
mdata        = mdata(3:end,2:end);
mnames       = mtxt(2,2:end);
fdata        = fdata(:,2:end);
end;

% Backcast macro series 46
y     = mdata(:,46);
dy    = y - mlags(y);
dX    = [ones(length(mdata),1),mdata(:,47:49)-mlags(mdata(:,47:49))];
ind   = find(isnan(dy)==0,1,'first');
reg   = nwest(dy(ind:end),dX(ind:end,:),0);
dy    = [dX(1:ind-1,:)*reg.beta;dy(ind:end)];
for t = ind-1:-1:1
   y(t) = y(t+1) - dy(t+1);
end

% Replace series 46 and apply Stock and Watson transformations
mdata(:,46) = y;
mdata       = mdata(13:end,:); %start in 1960
[yt,a,b,c]  = prepare(mdata,mnames,vartype);
yt          = yt(3:end,:); %remove NaNs
T           = length(yt);

% Save data to .mat file
data    = [yt,fdata(end-T+1:end,:)];
names   = [mnames,ftxt(1,2:end)]';
vartype = [vartype,repmat(9,1,size(fdata,2))];
dates   = 1900+(59:1/12:112-1/12)';
dates   = dates(end-T+1:end);
save jlndata data dates names vartype

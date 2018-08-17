function [yt,yh1,yh3,yh12]   = prepare(data,colheaders,vartype)
% -------------------------------------------------------------------------
% Transform raw data into stationary form
%   Input:  data       = raw data
%           colheaders = variable names
%           vartype    = variable type
%   Output: yt         = transformed data
% -------------------------------------------------------------------------

% Reformat data
series      = [];                                   %initialize variables
rawdata     = [];   
tcode       = [];
N           = size(data,2);                         %number of series
for i = 1:N                                         
    dum  = data(:,i);
    m    = mean(dum);
    if isnan(m) == 0;                               %check for non-numbers
        rawdata = [rawdata dum];                    %keep numerical data
        series  = str2mat(series,colheaders{1,i});  %store series name
        tcode   = [tcode; vartype(i)];              %store series type code
    else
        disp([i m]);                                %else show the error
    end
end
series   = series(2:end,:);                         %delete first name

% Transform data by group
group_id=[ones(2,1);ones(3,1)*4;ones(15,1); ...
          ones(30,1)*2;ones(10,1)*3;ones(10,1)*4;...
          ones(11,1)*5;ones(4,1)*8;ones(22,1)*6;...
          ones(21,1)*7; ones(3,1)*2; 4];
newtcode        = tcode;
newtcode(77)    = 7;  % compute percent change without taking logs
newtcode(51:60) = 5; % housing starts are differenced after logs
newtcode(46:49) = 2; % hours are differenced

y        = [];                                      %initialize y
yh1      = [];
yh3      = [];
yh12     = [];
N        = size(rawdata,2);                         %number of series kept
for i = 1:N
    % transxf.m is the new program that constructs the 
    % h=1,3,12 step ahead variables to be forecasted
    [dum,dum1,dum3,dum12] = transxf(rawdata(:,i),newtcode(i));
    y    = [y dum];
    yh1  = [yh1 dum1];
    yh3  = [yh3 dum3];
    yh12 = [yh12 dum12];
end
yt = y;

% Internal functions
 function [y,y1,y3,y12]=transxf(x,tcode)
 n=size(x,1);
 small=1e-6;
 y=NaN*ones(n,1);
 y1=NaN*ones(n,1);
 y3=NaN*ones(n,1);
 y12=NaN*ones(n,1);
 
switch(tcode);
  case 1,
  y=x;
  y1=lagn(x,-1);
  y3=lagn(x,-3);
  y12=lagn(x,-12);
 case 2,
  %y(1)=0;
  y(2:n)=x(2:n,1)-x(1:n-1,1);
  y1=(1200)*(lagn(x,-1)-x);
  y3=(1200/3)*(lagn(x,-3)-x);
  y12=(1200/12)*(lagn(x,-12)-x);
 case 3,
  %y(1)=0;y(2)=0;
  y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
  y1=1200*(lagn(x,-h)-x)-1200*(x-lagn(x,1));
  y3=(1200/3)*(lagn(x,-3)-x)-1200*(x-lagn(x,1));
  y12=(1200/12)*(lagn(x,-12)-x)-1200*(x-lagn(x,1));
 case 4,
  if min(x) < small; y=NaN; else;
  y=log(x);
  y1=lagn(log(x),-1);
  y3=lagn(log(x),-3);
  y12=lagn(log(x),-12);
  end;
 case 5,
  if min(x) > small; 
  x=log(x);
  %y(1)=0;
  y(2:n)=x(2:n)-x(1:n-1);
    y1=1200*(lagn(x,-1)-x);
  y3=(1200/3)*(lagn(x,-3)-x);
  y12=(1200/12)*(lagn(x,-12)-x);
  end;
 case 6,
  if min(x) > small;   
  %y(1)=0;y(2)=0;
  x=log(x);
  y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
  y1=1200*(lagn(x,-1)-x)-1200*(x-lagn(x,1));
  y3=(1200/3)*(lagn(x,-3)-x)-1200*(x-lagn(x,1));
  y12=(1200/12)*(lagn(x,-12)-x)-1200*(x-lagn(x,1));
  end;
case 7,
  %y1(1)=0;
  y1(2:n)=(x(2:n)-x(1:n-1))./x(1:n-1);
  %y(1)=0; y(2)=0;
  y(3:n)=y1(3:n)-y1(2:n-1);
  y1=1200*(lagn(x,-1)-x)./x-1200*(x-lagn(x,1))./lagn(x,1);
  y3=(1200/3)*(lagn(x,-3)-x)./x-1200*(x-lagn(x,1))./lagn(x,1);
  y12=(1200/12)*(lagn(x,-12)-x)./x-1200*(x-lagn(x,1))./lagn(x,1);
end;
 
  
y1(end-1+1,:) = NaN;
y3(end-3+1:end,:)=NaN;
y12(end-12+1:end,:)=NaN;
 end
end

function xx=trimr(x,a,b);
[nt,nc]=size(x);
if a >= 0 ;
xx=x(a+1:nt,:);
end;
if b >=0;
if a > 0; x=xx; end;
[nt,nc]=size(x);
xx=x(1:nt-b,1:nc);
end;
end

function y=lagn(x,n)
[nt,nc]=size(x);
if n> 0
x1=trimr(x,0,n);
 y=[zeros(n,nc); x1];
end
if n<0
x1=trimr(x,abs(n),0);
y=[x1;zeros(abs(n),nc)];
end
if n==0
  y=x;
end;
end 

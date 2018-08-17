function [x] = rshade(dates,p)
% -------------------------------------------------------------------------
% Shade NBER U.S. recession dates in a plot
%   Input:  dates = vector of sample dates (from 'datenum.m')
%           p     = shade bars on plot? 1 = yes, 0 = no;
%   Source: http://www.nber.org/cycles/cyclesmain.html
% -------------------------------------------------------------------------

% Options
if nargin==1;p=1;end;

% Input recession dates
R     = [datenum('15-Jun-1857'), datenum('15-Dec-1858');
         datenum('15-Oct-1860'), datenum('15-Jun-1861');
         datenum('15-Apr-1865'), datenum('15-Dec-1867');
         datenum('15-Jun-1869'), datenum('15-Dec-1870');
         datenum('15-Oct-1873'), datenum('15-Mar-1879');
         datenum('15-Mar-1882'), datenum('15-May-1885');
         datenum('15-Mar-1887'), datenum('15-Apr-1888');
         datenum('15-Jul-1890'), datenum('15-May-1891');
         datenum('15-Jan-1893'), datenum('15-Jun-1894');
         datenum('15-Dec-1895'), datenum('15-Jun-1897');
         datenum('15-Jun-1899'), datenum('15-Dec-1900');
         datenum('15-Sep-1902'), datenum('15-Aug-1904');
         datenum('15-May-1907'), datenum('15-Jun-1908');
         datenum('15-Jan-1910'), datenum('15-Jan-1912');
         datenum('15-Jan-1913'), datenum('15-Dec-1914');
         datenum('15-Aug-1918'), datenum('15-Mar-1919');
         datenum('15-Jan-1920'), datenum('15-Jul-1921');
         datenum('15-May-1923'), datenum('15-Jul-1924');
         datenum('15-Oct-1926'), datenum('15-Nov-1927');
         datenum('15-Aug-1929'), datenum('15-Mar-1933');
         datenum('15-May-1937'), datenum('15-Jun-1938');
         datenum('15-Feb-1945'), datenum('15-Oct-1945');
         datenum('15-Nov-1948'), datenum('15-Oct-1949');
         datenum('15-Jul-1953'), datenum('15-May-1954');
         datenum('15-Aug-1957'), datenum('15-Apr-1958');
         datenum('15-Apr-1960'), datenum('15-Feb-1961');
         datenum('15-Dec-1969'), datenum('15-Nov-1970');
         datenum('15-Nov-1973'), datenum('15-Mar-1975');
         datenum('15-Jan-1980'), datenum('15-Jul-1980');
         datenum('15-Jul-1981'), datenum('15-Nov-1982');
         datenum('15-Jul-1990'), datenum('15-Mar-1991');
         datenum('15-Mar-2001'), datenum('15-Nov-2001');
         datenum('15-Dec-2007'), datenum('15-Jun-2009')];
R   = year(R) + (month(R)-1)/12;
x   = zeros(size(dates));
for i = 1:length(R)
    x = x + ((R(i,1)<=dates) & (dates<=R(i,2)));
end

% Add bands to plot
if p ==1;
for r = 1:size(R,1);
    if R(r,2) >=min(dates)&&R(r,1)<=max(dates)
        ph(r) = patch([R(r,1) R(r,1) R(r,2) R(r,2)],...
                [get(gca,'Ylim') fliplr(get(gca,'Ylim'))],...
                [0,0,0,0],[.85 .85 .85]);
        box on;
        set(ph(r),'EdgeColor','none');
        uistack(ph(r),'bottom');
        set(gca,'layer','top');
        %h = findobj(gca,'Type','line');
        %set(h,'EraseMode','xor');
    end
end
end

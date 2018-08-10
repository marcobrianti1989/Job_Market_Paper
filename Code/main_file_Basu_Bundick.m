% Reading Data
filename                    = 'Stock_Prices';
sheet                       = 'Sheet1';
range                       = 'B1:E825';
data                        = xlsread(filename, sheet, range);
data                        = real(data);

j = 1;
for i = 1:length(data)
      if floor(((i - 1)/3)) == (i - 1)/3 && length(data) - i > 3
      SP(j,[1:4]) = mean(data(i:i+3,:),1);
      j = j + 1;
      end
end
plot(SP(:,1))
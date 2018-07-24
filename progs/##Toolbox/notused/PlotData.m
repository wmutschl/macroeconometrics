function PlotData(ENDO,VARopt)
figure;
for i = 1:VARopt.nvars    
    subplot(VARopt.nvars,1,i);    
    plot(datenum(VARopt.dates,'YYYYQQ'),ENDO(:,i));    
    datetick('x','YYYYQQ')
    title(VARopt.vnames(i))    
end
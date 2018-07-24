filename = 'C:\Users\w_muts01\Dropbox\shared folders\Statistik 1\daten\norwaydata.xlsm';
norwayfig = figure('Name','Data for Norway');
% GDP growth
gdpfig=figure();
[~,gdp_t,~] = xlsread(filename,'GDP','A5:A163');
gdp = xlsread(filename,'GDP','B5:B163');
cpi = xlsread(filename,'GDP','C5:C163');
gdpgrowth=((gdp(2:end)./cpi(2:end))./(gdp(1:end-1)./cpi(1:end-1))-1)*100;
%subplot(3,2,1)
plot(datenum(gdp_t(2:end),'QQ YYYY'),gdpgrowth,'LineWidth',1.5);
title('GDP growth (% - quarterly)','FontSize',26);
datetick('x','YYYY');
axis tight
xlim auto
set(gca,'fontsize',26);
set(gdpfig,'PaperOrientation','landscape');
set(gdpfig,'PaperUnits','normalized');
set(gdpfig,'PaperPosition', [0 0 1 1]);
print(gdpfig,'Norway_GDP_Growth','-dpdf','-r0');

% Unemployment Rate
urfig = figure();
[~,ur_t,~] = xlsread(filename,'Unemployment rate','A5:A155');
ur = xlsread(filename,'Unemployment rate','B5:B155');
%subplot(3,2,2)
plot(datenum(ur_t,'QQ YYYY'),ur,'LineWidth',1.5);  
datetick('x','YYYY');
title('Unemployment rate (rate - quarterly)','FontSize',26);
axis tight
xlim auto
set(gca,'fontsize',26);
set(urfig,'PaperOrientation','landscape');
set(urfig,'PaperUnits','normalized');
set(urfig,'PaperPosition', [0 0 1 1]);
print(urfig,'Norway_UR','-dpdf','-r0');

% 3-month interest rate
intfig = figure();
[~,int_t,~] = xlsread(filename,'3 month interest rates','A5:A328');
int = xlsread(filename,'3 month interest rates','E5:E328');
%subplot(3,2,3)
plot(datenum(int_t,'dd.mm.yyyy'),int,'LineWidth',1.5);  
datetick('x','YYYY');
title('3-month interest rate (rate - monthly)','FontSize',26);
axis tight
xlim auto
set(gca,'fontsize',26);
set(intfig,'PaperOrientation','landscape');
set(intfig,'PaperUnits','normalized');
set(intfig,'PaperPosition', [0 0 1 1]);
print(intfig,'Norway_int','-dpdf','-r0');

% Oslo Stock Exchange
stockfig = figure();
[~,stock_t,~] = xlsread(filename,'Stock Exchange','A5:A9150');
stock = xlsread(filename,'Stock Exchange','B5:B9150');
%subplot(3,2,4)
plot(datenum(stock_t,'dd.mm.yyyy'),stock,'LineWidth',1.5);  
title('Oslo Stock Exchange (index - daily)','FontSize',26);
datetick('x','YYYY');
axis tight
xlim auto
set(gca,'fontsize',26);
set(stockfig,'PaperOrientation','landscape');
set(stockfig,'PaperUnits','normalized');
set(stockfig,'PaperPosition', [0 0 1 1]);
print(stockfig,'Norway_stock','-dpdf','-r0');

% Population
popfig = figure();
[pop_t,~,~] = xlsread(filename,'Population','A5:A41');
pop = xlsread(filename,'Population','B5:B41');
%subplot(3,2,5)
plot(pop_t,pop./1000000,'LineWidth',1.5);  
title('Population (millions - yearly)','FontSize',26);
axis tight
xlim auto
set(gca,'fontsize',26);
set(popfig,'PaperOrientation','landscape');
set(popfig,'PaperUnits','normalized');
set(popfig,'PaperPosition', [0 0 1 1]);
print(popfig,'Norway_pop','-dpdf','-r0');

% Real house prices
housefig = figure();
[~,house_t,~] = xlsread(filename,'Residential Property Prices','A37:A195');
house = xlsread(filename,'Residential Property Prices','B37:B195');
%subplot(3,2,6)
plot(datenum(house_t,'QQ YYYY'),house./cpi,'LineWidth',1.5);  
title('Real house prices (index - quarterly)','FontSize',26);
datetick('x','YYYY');
axis tight
xlim auto
set(gca,'fontsize',26);
set(housefig,'PaperOrientation','landscape');
set(housefig,'PaperUnits','normalized');
set(housefig,'PaperPosition', [0 0 1 1]);
print(housefig,'Norway_house','-dpdf','-r0');

% Export
figure();
suptitle('Data For Norway')
set(gcf, 'Color', 'w');
print(norwayfig,'norwaydata','-dpdf','-r0','-fillpage')
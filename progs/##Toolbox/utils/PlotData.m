function PlotData(ENDO,EXOG,opt)
nvars = size(ENDO,1);
nvars_ex = size(EXOG,1);
dates_freq = opt.dates_freq;
dates = opt.dates;
vnames = opt.vnames;
vnames_ex = opt.vnames_ex;

% Plot dataset for endogenous variables
figure;
for i = 1:nvars
    subplot(nvars,1,i);
    if strcmp(dates_freq,'m')
        plot(datenum(dates,'YYYYmm'),ENDO(:,i));    
        datetick('x','YYYYmm');
    elseif strcmp(dates_freq,'Q')
        plot(datenum(dates,'YYYYQQ'),ENDO(:,i));    
        datetick('x','YYYYQQ');
    end
    title(vnames(i));
end

% Plot dataset for exogenous variables
if nvar_ex > 0
    figure;
    for i = 1:nvars
        subplot(nvars,1,i);
        if strcmp(dates_freq,'m')
            plot(datenum(dates,'YYYYmm'),ENDO(:,i));    
            datetick('x','YYYYmm');
        elseif strcmp(dates_freq,'Q')
            plot(datenum(dates,'YYYYQQ'),ENDO(:,i));    
            datetick('x','YYYYQQ');
        end
        title(vnames(i));
    end
end
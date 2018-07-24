function FEVDplot(opt,FEVDpoint,FEVDINF,FEVDSUP)
% =======================================================================
% Plot the FEVDs computed with VARfevd
% =======================================================================
% VARfevdplot(FEVD,VARopt,vnames,INF,SUP)
% -----------------------------------------------------------------------
% INPUT
%   - FEVD(:,:,:): matrix with 't' steps, the FEVD due to 'j' shock for 
%       'k' variable
%	- VARopt: options of the FEVDs (see VARoption)
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%   - INF: lower error band
%   - SUP: upper error band
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com


%% Check inputs
%================================================
vnames = opt.vnames;

%% Retrieve and initialize variables 
%================================================
filename = ['figures/', opt.dataset '_FEVD_'];

% Initialize FEVD matrix
[nsteps, nvars, nshocks] = size(FEVDpoint);

% Define the rows and columns for the subplots
row = round(sqrt(nshocks));
col = ceil(sqrt(nshocks));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);



%% Plot
%================================================
figure('units','normalized','outerposition',[0 0.1 1 0.9])
for ii=1:nvars
    for jj=1:nshocks
        subplot(row,col,jj);
        plot(steps,FEVDpoint(:,jj,ii),'LineStyle','-','Color',[0.01 0.09 0.44],'LineWidth',2);
        hold on
        plot(x_axis,'k','LineWidth',0.5)
        if isempty(FEVDINF)==0 && isempty(FEVDSUP)==0
            plot(steps,FEVDINF(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
            hold on
            plot(steps,FEVDSUP(:,jj,ii),'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
        end
        xlim([1 nsteps]); ylim([0 1]);
        title([vnames{ii} ' to ' vnames{jj}], 'FontWeight','bold','FontSize',10); 
    end
    % Save
    FigName = [filename num2str(ii)];
    print('-dpng','-r100',FigName);
    print('-dpdf','-r100',FigName);
    clf('reset');
end

close all

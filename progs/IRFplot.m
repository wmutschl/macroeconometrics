function IRFplot(opt,IRFpoint,IRFINF,IRFSUP)
% =======================================================================
% Plot the IRFs computed with VARir
% =======================================================================
% VARirplot(IRF,VARopt,vnames,INF,SUP)
% -----------------------------------------------------------------------
% INPUT
%   - IRF(:,:,:) : matrix with periods, variable, shock
%   - VARopt: options of the VAR (see VARopt from VARmodel)
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
IRFcumsum = opt.IRFcumsum;
%% Retrieve and initialize variables 
%================================================
filename = ['./figures/', opt.dataset, '_IRF_'];

% Initialize IRF matrix
[nsteps, nvars, nshocks] = size(IRFpoint);

% Define the rows and columns for the subplots
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);


%% Plot
%================================================
figure('units','normalized','outerposition',[0 0.1 1 0.9]);
for jj=1:nshocks                
    for ii=1:nvars
        subplot(row,col,ii);
        irfpoint = IRFpoint(:,ii,jj);
        if IRFcumsum 
            irfpoint = cumsum(irfpoint);
        end
        plot(steps,irfpoint,'LineStyle','-','Color',[0.01 0.09 0.44],'LineWidth',2);
        hold on
        plot(x_axis,'k','LineWidth',0.5)
        if isempty(IRFINF) == 0 && isempty(IRFSUP)==0
            irfinf = IRFINF(:,ii,jj);
            irfsup = IRFSUP(:,ii,jj);
            if IRFcumsum
                irfinf = cumsum(irfinf);
                irfsup = cumsum(irfsup);
            end
            plot(steps,irfinf,'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
            hold on
            plot(steps,irfsup,'LineStyle',':','Color',[0.39 0.58 0.93],'LineWidth',1.5);
        end
        xlim([1 nsteps]);
        title([vnames{ii} ' to ' vnames{jj}], 'FontWeight','bold','FontSize',10); 
    end
    % Save
    FigName = [filename num2str(jj)];
    print('-dpng','-r100',FigName);
    print('-dpdf','-r100',FigName);
    clf('reset');
end

close all

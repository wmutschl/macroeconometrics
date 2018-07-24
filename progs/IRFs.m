function IRFpoint  = IRFs(Acomp,impact,opt)
% =======================================================================
% Compute IRFs for a SVAR model given the companion matrix of the
% reduced-form VAR and the impact matrix
% =======================================================================
% IRFpoint  = IRFs(Acomp,B0inv,p,H,opt)
% -----------------------------------------------------------------------
% INPUT
%   - Acomp: Companion matrix of reduced-form VAR
%   - impactMat: Impact matrix of SVAR model, either inv(B0) or inv(B0)*sqrt(E[epsilon_t*epsilon_t'])
%   - opt: structure of options
%   -      - vnames: cell of strings with names for the variables to plot
%          - epsnames: cell of strings with names for the structural shocks
%          - IRFcumsum: use cumulative sum for IRFs (1) or not (0), vector nvars x 1
%          - doplot: 1 displays the plots
%          - dosave: 1 saves the plots as png or pdf
%          - filename: name of file to save plots as png and pdf
% ----------------------------------------------------------------------- 
% OUTPUT
%   - IRFpoint(j,k,h): matrix with h=0,1,...,nsteps+1 steps, containing the IRF of 'j' variable to 'k' shock
% =======================================================================
% Willi Mutschler, November 2017
% willi@mutschler.eu

% Get options
nlag   = opt.nlag;
nsteps = opt.nsteps;
nvars  =size(impact,1);
% Set default options if not specified in opt structure
try vnames = opt.vnames; catch; vnames = sprintfc('y%d',1:nvars);end
try epsnames = opt.epsnames; catch; epsnames = sprintfc('eps%d',1:nvars);end
try IRFcumsum = opt.IRFcumsum; catch; IRFcumsum=zeros(1,nvars); end
try doplot = opt.doplot; catch; doplot=1; end
try dosave = opt.dosave; catch; dosave=0; end
try filename = [opt.filename '_IRF_']; catch; filename = 'IRF_'; end

% Initialize variables 
IRFpoint = nan(nvars,nvars,nsteps+1);
J = [eye(nvars) zeros(nvars,nvars*(nlag-1))];
JtB0inv = J'*impact;

% Compute the impulse response
AA = eye(size(Acomp));
for h=1:(nsteps+1)    
    IRFpoint(:,:,h) = J*AA*JtB0inv;
    AA = AA*Acomp;
end

% If we need to cumsum
for ivar = 1:nvars
    if IRFcumsum(ivar) == 1
        IRFpoint(ivar,:,:) = cumsum(IRFpoint(ivar,:,:),3);
    end
end

% Plot
if doplot
    % Define a timeline
    steps = 0:1:nsteps;
    x_axis = zeros(1,nsteps+1);
    figure('units','normalized','outerposition',[0 0.1 1 0.9]);
    count = 1;
    for ishocks=1:nvars % index for
        for ivars=1:nvars % index for
            irfpoint = squeeze(IRFpoint(ivars,ishocks,:));
            subplot(nvars,nvars,count);
            plot(steps,irfpoint,steps,x_axis,'k','LineWidth',2);
            xlim([0 nsteps]);
            title([epsnames{ishocks}, ' Shock'], 'FontWeight','bold','FontSize',10);
            ylabel(vnames{ivars}, 'FontWeight','bold','FontSize',10);
            count = count+1;
        end
    end
    if dosave
        % Save
        FigName = [filename num2str(ivars)];
        print('-dpng','-r0',FigName);
        print('-dpdf','-r0','-bestfit',FigName);
    end
end
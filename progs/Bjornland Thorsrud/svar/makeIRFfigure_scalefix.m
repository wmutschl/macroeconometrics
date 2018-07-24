function makeIRFfigure_scalefix(irf_m,irf_q,name,id,legendNames)
% PURPOSE: Make impulse response figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nirf=size(irf_m,1);

f=figure('name',['IRF - r' name{1} '-s' name{2} '-' id]);
h=plot(irf_m,'k','linewidth',2);
hold on
h1=plot(irf_q,'k--','linewidth',2);
h2=plot(zeros(nirf,1),'k','linewidth',1);

set(gca,'xlim',[1 nirf]);
%set(gca,'xticklabel',num2str(str2num(get(gca,'xticklabel'))-1,'%6.0f'))
ylim=get(gca,'ylim');
set(gca,'ylim',[ylim(1) ylim(2)+1]);

if nargin==5
    lh=legend([h h1],{legendNames{1},legendNames{2}});
    set(lh,'box','off','location','northwest');    
end;



function ps=printSettings(h,hndl)
%% PURPOSE: Define default print settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ps.lineWidth=2;
ps.fontSize=12;
ps.figureFormat.x={'eps'};%pdf
ps.transparentFigures=true;
ps.renderer='printer'; %'opengl'

if nargin>=1
    set(h,'fontsize',ps.fontSize);
    set(h,'box','off');              
end;
if nargin==2
    scrsz=get(hndl,'Position');
    set(hndl,'Position',[scrsz(1) scrsz(2) scrsz(3)*2 scrsz(4)]);            
end;




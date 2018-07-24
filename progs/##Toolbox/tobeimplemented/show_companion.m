function show_companion(beta,alpha,Psi,tm,curr_ini,ind,rind)

file_str = '';
res_str = '';
if ind==1;
   file_str = 'rest-';
   res_str = ' restricted';
end;
realscrsz = get(0,'ScreenSize');
%
if rind==1;
   eigenvalues = companion(tm.Z0,tm.Z1,tm.Z2,tm.r,tm.p,tm.n,tm.q1,tm.q0,tm.d,tm.T);
else;
   eigenvalues = companionrest(alpha,Psi,beta,tm.p,tm.n,tm.d);
end;

eigenvalues = sort(eigenvalues);
eigenvalues = [real(eigenvalues) imag(eigenvalues)];
if curr_ini.colors==1;
   colorstr = 'r';
   colorstr2 = 'b';
else;
   colorstr = [0.5 0.5 0.5];
   colorstr2 = 'k';
end;
myplot = figure('Units','pixels', ...
   'Position',[(realscrsz(3)-550)/2 (realscrsz(4)-360)/2 550 360], ...
   'Menubar','none', ...
   'Toolbar','none', ...
   'NumberTitle','off', ...
   'Visible','off', ...
   'Name','Eigenvalues of the Companion Matrix');
%
% obey centering and filling defaults
%
setfigureproperties(myplot);
%
%
% Plot unit circle
%
my_x = pi*[0:.5:2];
my_y = [0 1 0 -1 0 1 0; 1 0 1 0 -1 0 1];
my_pp = spline(my_x,my_y);
yy = ppval(my_pp, linspace(0,2*pi,101));
plot(yy(1,:),yy(2,:),'-b','Color',colorstr2,'Marker','none','LineWidth',1,'UserData','noshow');
hold('on');
plot(zeros(2,1),[min([min(eigenvalues(:,2)) -1])-0.1 max([max(eigenvalues(:,2)) 1])+0.1]','LineWidth',1,'Color','k');
hold('on');
plot([min([min(eigenvalues(:,1)) -1])-0.1 max([max(eigenvalues(:,1)) 1])+0.1]',zeros(2,1),'LineWidth',1,'Color','k');
hold('on');
plot(eigenvalues(:,1),eigenvalues(:,2),'LineStyle','none','Marker','+','LineWidth',1,'Color',colorstr);
%
% Make axes into a square box
%
set(gca,'Units','pixels');
thispos = get(gca,'Position');
if iscompiled==1;
   addnum = 22;
else;
   addnum = 0;
end;
newleft = (550-thispos(4))/2+(addnum/2);
set(gca,'Position',[newleft thispos(2) thispos(4)-addnum thispos(4)]);
set(gca,'Units','normalized');
set(gca,'XLim',[min([min(eigenvalues(:,1)) -1])-0.1   max([max(eigenvalues(:,1)) 1])+0.1]);
set(gca,'YLim',[min([min(eigenvalues(:,2)) -1])-0.1   max([max(eigenvalues(:,2)) 1])+0.1]);
setaxesfonts('Title',['Eigenvalues of companion matrix for' res_str ' rank ' num2str(tm.r,'%0.0f')],'Xlabel','Real','YLabel','Imaginary');
toolbar_fcn(gca);
hold('off');
if curr_ini.external==0;
   set(myplot,'Visible','on');
else;
   file = [tm.outputdir 'eigenvalues-' file_str num2str(tm.r)];
   mstat = autosavefigure(myplot,file);
   if mstat==0;
      autoshowfigure(file);
   end;
end;

%
% end of show_companion.m
%

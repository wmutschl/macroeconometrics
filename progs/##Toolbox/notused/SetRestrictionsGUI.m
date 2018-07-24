function SetRestrictionsGui(opt)
nvars = opt.nvars;
short = opt.RESTR.short;
long = opt.RESTR.long;
sign = opt.RESTR.sign;

close all
arr_fig = figure('unit','normalized','NumberTitle','off','Menubar','none','resize','on','Name','Restrictions');
ok_but = uicontrol('Style','pushbutton','unit','normalized','String','OK','position',[0.45 0 .1 .1],'callback','uiresume','tag','ok');

posx = .05;
posxb = posx;
posy = 1-.15;

uicontrol('Style','text','unit','normalized','position', [0 0.5 1 .5],'String','Please set Restrictions on inv(B0)');
for i = 1:nvars
    for j = 1:nvars
        a(i,j)=uicontrol('Style','edit','unit','normalized','position', [posx posy .05 .05],'String',short(i,j));
        posx = posx+0.05;
    end
    posx = posxb;
    posy = posy-0.05;
end


uicontrol('Style','text','unit','normalized','position', [0 0.5 1 .5],'String','Please set Restrictions on inv(B0)');
for i = 1:nvars
    for j = 1:nvars
        a(i,j)=uicontrol('Style','edit','unit','normalized','position', [posx posy .05 .05],'String',short(i,j));
        posx = posx+0.05;
    end
    posx = posxb;
    posy = posy-0.05;
end


movegui('center');
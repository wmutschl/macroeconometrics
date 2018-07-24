function data = mltable(fig, hObj, action, columnInfo, rowHeight, cell_data, gFont, varargin)
% function data = mltable(fig, hObj, action, columnInfo, rowHeight, cell_data, gFont)
%
% Author: Morris Maynard
% Based on code by Gregory Gershanok
% Last update: 27 Jan 2005
%
% Manages a table with editing and scrolling capability
% Features: varying column widths, entry formatting, insert/delete rows,
%           editable or read-only cells, font control, scaled numeric
%           display, row selection highlighing, multiple tables per figure,
%           optional checkboxes on left-hand side
%
% Usage:
% Supply the parent figure, the handle of the axes object to use, the
% 'CreateTable' action, info about columns, and the cell data 
%
% Example usage: (also just run with no arguments to see result)
% 
% fig = nf;
% tbl = axes('units', 'pixels','position', [10 10 400 100]);
% cell_data = {... 
%           'Alpha',   1, 2, 3,'';...
%           'Bravo',   4, 5, 6,'';...
%           'Charlie', 7, 8, 9,'';...
%           'Dog',    10,11,12,'';...
%           'Echo',   13,14,15,'';...
%           'Foxtrot',16,17,18,'';...
%           'Golf',   19,20,21,'';...
%           'Hotel',  26,27,28,'';...
%           };
% 
% columninfo.titles={'Param','Lower Limit','Upper Limit','Initial Value','Result'};
% columninfo.formats = {'%4.6g','%4.6g','%4.6g','%4.6g', '%4.6g'};
% columninfo.weight =      [ 1, 1, 1, 1, 1];
% columninfo.multipliers = [ 1, 1, 1, 1, 1];
% columninfo.isEditable =  [ 1, 1, 1, 1, 0];
% columninfo.isNumeric =   [ 0  1, 1, 1, 1];
% columninfo.withCheck = true; % optional to put checkboxes along left side
% columninfo.chkLabel = 'Use'; % optional col header for checkboxes
% rowHeight = 16;
% gFont.size=9;
% gFont.name='Helvetica';
% 
% mltable(fig, tbl, 'CreateTable', columninfo, rowHeight, cell_data, gFont);
%
% To use in a GUIDE-created figure:
%
% Create a figure including a blank "axes" object with the tag 'tblParams'
% Put the lines starting with the "cell_data" line above into your figure's
% OpeningFcn, but replace the mltable line with:
%
% mltable(gcf, handles.tblParams, 'CreateTable', columninfo, rowHeight, cell_data, gFont);
%% so clicking outside the table will finish edit in progress...
% endfcn = sprintf('mltable(%14.13f, %14.13f, ''SetCellValue'');', hObject, handles.tblParams);
% set(hObject,'buttondownfcn',endfcn);
%
% To access the data edited by the table:
%
% info = get(tbl, 'userdata');
% data = info.data;
%

%-------------------------------------------------------------------------
% All functions dispatched from here. 
% If necessary to call from figure use: mltable(fig, hObj, 'Action',...)

global MINROWS;
MINROWS = 3;
if ~exist('action', 'var')
    data = mltest;
    return
end
    
switch(action)
  case 'CreateTable'
    data = createTable(fig, hObj, columnInfo, rowHeight, cell_data, gFont);
  case 'DestroyTable'
    data = destroyTable(fig, hObj);
  case 'ResizeTable'
    fig = resizeTable(fig, hObj);
  case 'ScrollData'
    fig = scrollData(fig, hObj);
  case 'EditCell'
    editCell(fig, hObj);
  case 'SetCellValue'
    setCellValue(fig, hObj);
  case 'EndEdit'
    setCellValue(fig, hObj);
  case 'AddRow'
    addRow(fig, hObj);
  case 'DelRow'
    delRow(fig, hObj);
  case 'SetOnSetCell'
    setOnSetCell(fig, hObj, varargin{:});
  case 'OnSetCell'
    data = onSetCell(fig, hObj, varargin{:});
  case 'SetDblClick'
    setDblClick(fig, hObj, varargin{:});
  case 'OnDblClick'
    onDblClick(fig, hObj, varargin{:});
  case 'OnCheck'
    onCheck(fig, hObj, varargin{1});
  case 'SetCheck'
    setCheck(fig, hObj, varargin{1}, varargin{2});
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = createTable(fig, hObj, columnInfo, rowHeight, cell_data, gFont)
% Initially creates the table
%-------------------------------------------------------------------------
global MINROWS
data.figure = fig;
% get axes position in pixel coordinates
set(hObj, 'units', 'pixels');
set(hObj, 'visible', 'on');
pos_ax = get(hObj, 'position');
% set up grid info structure
ds = size(cell_data);
data.maxRows = ds(1);
if data.maxRows < MINROWS
    blanks = cell(1, ds(2));
    for ii = data.maxRows+1:MINROWS
        cell_data = [cell_data; blanks];
    end
    data.maxRows = MINROWS;
end
data.data = cell_data;
data.isChecked = zeros(1,size(cell_data, 1));
data.axes = hObj;
data.userModified = zeros(ds);
data.rowHeight = rowHeight;
data.columnInfo = columnInfo;
data.numCols= length(columnInfo.titles);
data.ltGray = [92 92 92]/255;
data.OffscreenPos = [-1000 -1000 30 20];
data.selectedRow = 0;
data.selectedCol = 0;
data.gFont = gFont;

data.doCheck = false;
if isfield(data.columnInfo, 'withCheck') && ...
    data.columnInfo.withCheck ~= 0
    data.doCheck = true;
end

% use 0...1 scaling on table x and y positions
set(fig, 'CurrentAxes', data.axes);
set(data.axes, 'box', 'on', 'DrawMode', 'fast');
set(data.axes, 'xlimmode', 'manual', 'xlim', [0 1], 'ylim', [0 1], ...
               'xtick', [], 'ytick', [], 'xticklabelmode', 'manual', 'xticklabel', []);
           
if data.doCheck % shrink on left for checkboxes column
    data.checkdx = pos_ax(3) * 20 * (1/pos_ax(3)); % chkbox offset
    pos_ax(1) = pos_ax(1) + data.checkdx;
    pos_ax(3) = pos_ax(3) - data.checkdx;
end
pos_ax(3) = pos_ax(3) - 10; % width of slider
set(data.axes, 'position', pos_ax, 'LineWidth', 2);
% callback for starting editing 
editfcn = sprintf('mltable(%14.13f, %14.13f, ''EditCell'');',fig, hObj);
set(data.axes, 'ButtonDownFcn', editfcn);
% callback for scrolling table
scrfcn = sprintf('mltable(%14.13f, %14.13f, ''ScrollData'');',fig, hObj);
data.slider = uicontrol('style', 'slider', 'units', 'pixels',...
    'position', [pos_ax(1)+pos_ax(3)+2 pos_ax(2) 16 pos_ax(4)],...
    'Callback', scrfcn);

% Add buttons for addrow/delrow
if sum(columnInfo.isEditable) > 0 && (~isfield(columnInfo,'rowsFixed') ||...
        ~columnInfo.rowsFixed)
	btnw = 19; btnh = 15;
	btnx = pos_ax(1) + pos_ax(3) - btnw - 2;
	btny = pos_ax(2) + pos_ax(4) + 2;
	btnfcn = sprintf('mltable(%14.13f, %14.13f, ''AddRow'');',fig, hObj);
	data.btnAdd = uicontrol('style', 'pushbutton', 'units', 'pixels',...
        'position', [btnx, btny, btnw, btnh],...
        'string',' + ','fontsize',12,'Callback', btnfcn,...
        'TooltipString','Click to add a row');
	btnfcn = sprintf('mltable(%14.13f, %14.13f, ''DelRow'');',fig, hObj);
	data.btnDel = uicontrol('style', 'pushbutton', 'units', 'pixels',...
        'position', [btnx + btnw + 2, btny, btnw, btnh],...
        'string',' - ','fontsize',12,'Callback', btnfcn,...
        'TooltipString','Click to remove selected row');
	
	set(data.btnAdd,'Units','normalized');
	set(data.btnDel,'Units','normalized');
else
    data.btnAdd = [];
    data.btnDel = [];
end

set(hObj, 'UserData', data);
% so clicking outside the table will finish edit in progress
endfcn = sprintf('mltable(%14.13f, %14.13f, ''SetCellValue'');', fig, hObj);
set(fig,'buttondownfcn',endfcn);

resizeTable(fig, hObj);

return

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = resizeTable(fig, hObj)
% fit table within boundaries and update scrollbar
% at this time doesn't handle figure resize
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if isempty(data)
    return
end
% zap checkboxes if any
if data.doCheck && isfield(data, 'checks')
    delete(data.checks(:));
    data.checks = [];
end

cla(hObj);

set(hObj, 'units', 'pixels');
set(fig,'CurrentAxes',hObj);  

pos_ax = get(hObj,'position');
data.numRows = floor((pos_ax(4)-(2*data.rowHeight))/data.rowHeight);
if data.numRows > data.maxRows
    data.numRows = data.maxRows;
end

unit_d_h = 1/pos_ax(3);
unit_d_v = 1/(pos_ax(4) - data.rowHeight);

if(data.numRows < data.maxRows)
    set(data.slider,'Units','pixels');
	set(data.slider, 'visible', 'on',...
                   'position', [pos_ax(1)+pos_ax(3)+1 pos_ax(2) 16 pos_ax(4)], ...
                   'min', 0,...
                   'max', data.maxRows - data.numRows, 'value', data.maxRows - data.numRows, ...
                   'sliderstep', [1/(data.maxRows-data.numRows) data.numRows/(data.maxRows-data.numRows)]);
else  
  set(data.slider, 'visible', 'off');
end

% get gui units for rows and columns
% average column width and row height
d_x = 1/sum(data.columnInfo.weight);
d_y = 1/(data.numRows + 1);
% minimum adjust unit
unit_d_h = 1/(pos_ax(3));
unit_d_v = 1/(pos_ax(4) - data.rowHeight);

% Horizontal line positions
lx_h = ones(2, data.numRows+1);
lx_h(1, :) = 0;
ly_h = [1:data.numRows+1; 1:data.numRows+1]/(data.numRows+1);

% Vertical line positions
ly_v = ones(2, data.numCols);
ly_v(1, :) = 0;
lx_v = [d_x*data.numCols 2:data.numCols; d_x*data.numCols 2:data.numCols]/data.numCols;
for i = 2:data.numCols
  lx_v(:, i) = lx_v(:, i-1)+d_x*data.columnInfo.weight(i);
end
% draw initial grid
data.vertLines  = line(lx_v, ly_v);
data.vertLines1  = line(lx_v(:, 1:(data.numCols-1)), ly_v(:, 1:(data.numCols-1)));
data.horizLines = line(lx_h, ly_h);
set(data.horizLines, 'color', data.ltGray);
set(data.vertLines, 'color', data.ltGray, 'LineWidth', 2);
set(data.vertLines1, 'color', [1 1 1], 'LineWidth', 0.5);

% now display text in grid     
txt_x = [0:data.numCols-1]/data.numCols + 4*unit_d_h;
for i = 2:data.numCols
  txt_x(i) = txt_x(i-1) + d_x*data.columnInfo.weight(i-1);
end
data.txt_x = txt_x;
data.txtCells = zeros(data.numRows, data.numCols);
uictx = get(hObj,'UIContextMenu');

chkdx = (pos_ax(3) * 20 * unit_d_h); % chkbox offset
chkdy = pos_ax(4) * d_y;

for j = 1:data.numRows
	txt_y = (data.numRows-j+1)/(data.numRows+1) * ones(1, data.numCols);
	txt_y = txt_y - (0.9*d_y); % reduce by (almost?) one row (title row)
	data.txtCells(j, :) = text(txt_x, txt_y, 'a','Clipping','on');
	if data.doCheck % put checkboxes to left of row
		poschk = [ pos_ax(1) - chkdx...
                 pos_ax(2) + 3 + chkdy * (data.numRows - j)...
                 10 10 ];
		data.checks(j) = ...
          uicontrol('style','checkbox','units','pixels','position', poschk);
		chkfcn = sprintf('mltable(%14.13f, %14.13f,''OnCheck'',[], [], [], [], %d);',...
          fig, hObj, j);
		set(data.checks(j), 'Units','normalized','Callback', chkfcn);
	end
	for i = 1:data.numCols
		if data.columnInfo.isNumeric(i)
			nums = data.data{j, i}/data.columnInfo.multipliers(i);
			nums = num2str(nums, data.columnInfo.formats{i});
			set(data.txtCells(j, i), 'string', nums);
		else
            set(data.txtCells(j, i), 'string', StripChars(data.data{j, i}));
		end
		set(data.txtCells(j, i), 'Position', [txt_x(i), txt_y(i)]);
		set(data.txtCells(j, i), 'UIContextMenu', uictx);
	end
end
set(data.txtCells(:, :), 'FontSize', data.gFont.size, 'FontName', data.gFont.name, 'FontWeight', 'normal', ...
                         'HorizontalAlignment','left','VerticalAlignment', 'bottom');
set(data.txtCells(1:data.numRows, 1:data.numCols), 'buttondownfcn', get(data.axes, 'ButtonDownFcn'));

% do title cells
titleCellsy = ones(1, data.numCols);
data.titleCells = text(txt_x, titleCellsy, data.columnInfo.titles);
set(data.titleCells(:, :), 'FontSize', data.gFont.size, 'FontName', data.gFont.name, 'FontWeight', 'bold', ...
                         'HorizontalAlignment','Center','VerticalAlignment', 'top','Editing','off');
for i = 1:length(data.titleCells)
	if i > 1
        titleOffset = 0.92 * (lx_v(1,i) - lx_v(1,i-1)) / 2;
	else
        titleOffset = (lx_v(1,i)) / 2;
	end
	pos = get(data.titleCells(i), 'position');
	pos(1) = pos(1) + titleOffset;
	set(data.titleCells(i), 'position', pos);
	if(data.columnInfo.isEditable(i) == 0)
		set(data.titleCells(i, :), 'FontWeight', 'normal');
		set(data.txtCells(:, i), 'FontWeight', 'normal');
	end
end
if data.doCheck && isfield(data.columnInfo,'chkLabel')
    chkLbl = text(-25*unit_d_h, titleCellsy(1), data.columnInfo.chkLabel);
    set(chkLbl, 'FontSize', data.gFont.size, 'FontName', data.gFont.name, 'FontWeight', 'bold', ...
               'HorizontalAlignment','left','VerticalAlignment', 'top','Editing','off');
end

row = 1;
x1 = 0; x2 = 1;
y1 = (data.numRows - row + 1)/(data.numRows+1);
y2 = y1 + (1/(data.numRows+1)) + .01;
makepatch(hObj,x1,x2,y1,y2, 0.753);

set(hObj, 'UserData', data);

scrollData(fig, hObj);
set(hObj, 'units', 'normalized');
return

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = scrollData(fig, hObj)
% handle scrollbar (slider) callback
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if isempty(data)
    return
end

set(data.slider,'Units','pixels');
if isfield(data,'editBox') && ishandle(data.editBox)
    delete(data.editBox);
end

% handle non-scroll case in case slider was switched off
if(strcmp(get(data.slider, 'visible'), 'off') == 1)
	ind0 = 0;
	for i = 1:data.numRows
        if data.doCheck
            set(data.checks(i),'value',data.isChecked(i+ind0));
        end
		for j = 1:data.numCols
            if data.columnInfo.isNumeric(j)
                nums = data.data{i, j}/data.columnInfo.multipliers(j);
                nums = num2str(nums, data.columnInfo.formats{j});
                set(data.txtCells(i, j), 'string', nums);
            else
                set(data.txtCells(i, j), 'string', data.data{i, j});
            end
		end
	end
else	
	val = get(data.slider, 'Value');
	max_val = get(data.slider, 'Max');
	min_val = get(data.slider, 'Min');
	
	val0 = data.maxRows - data.numRows;
	ind0 = round(val0-val);
	% move the text to give illusion of scrolling
	for i = ind0+1:(ind0+data.numRows)
      if data.doCheck
          set(data.checks(i-ind0),'value',data.isChecked(i));
      end
      for j = 1:data.numCols
            if data.columnInfo.isNumeric(j)
                nums = data.data{i, j}/data.columnInfo.multipliers(j);
                nums = num2str(nums, data.columnInfo.formats{j});
                set(data.txtCells(i-ind0, j), 'string', nums);
            else
                set(data.txtCells(i-ind0, j), 'string', data.data{i, j});
            end
      end
	end
end
% save scroll position
data.ind0 = ind0;

data.hpatch = remakepatch(data.selectedRow, data, hObj);

set(data.slider,'Units','normalized');
set(hObj, 'UserData', data);
return

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function [hpatch] = makepatch(hObj, x1, x2, y1, y2, co)
% -----------------------------------------------------------------------
    if ~exist('co','var')
        co = 0.88;
    end
    hpatch = ...
        patch([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],[co co co],...
        'HitTest','off','FaceAlpha',0.2);
    ch=get(hObj,'children');
    % put new patch at bottom layer
    c1 = ch(1); % latest object
    ch = ch(2:end); % shift up
    ch(size(ch,1)+1) = c1; % append new
   set(hObj,'children',ch);
return

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function [hpatch] = remakepatch(row, data, hObj)
% -----------------------------------------------------------------------
hpatch = []; % empty return if nothing to do
% see if selected row is visible
if ((row - data.ind0) <= (data.numRows) && (row - data.ind0) > 0)
    % yes, compute row coordinates
    row = row - data.ind0 + 1;
    x1 = 0; x2 = 1;
    y1 = (data.numRows - row + 1)/(data.numRows+1);
    y2 = y1 + (1/(data.numRows+1)) + .005;
    % see if a previous patch exists
	if isfield(data,'hpatch') && ~isempty(data.hpatch) && ishandle(data.hpatch)
        % yes, see if it is on same row
        yp = get(data.hpatch,'ydata');
        if y1 ~= yp(1)
            % no, delete old patch and make new one
            delete(data.hpatch);
            data.hpatch = makepatch(hObj,x1,x2,y1,y2);
        end
	else % no previous patch exists, make new one
		data.hpatch = makepatch(hObj,x1,x2,y1,y2);
	end
    hpatch = data.hpatch; % return new or previous patch
else % if patch is no longer visible, delete it
	if isfield(data,'hpatch') && ~isempty(data.hpatch) && ishandle(data.hpatch)
        delete(data.hpatch);
    end
end
return

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = destroyTable(fig, hObj)
% Destroys the table
%-------------------------------------------------------------------------
    data = get(hObj, 'UserData');
    if ~isempty(data)
        set(hObj, 'visible','off');
        cla(hObj);
        
		if isfield(data.columnInfo, 'withCheck') && ...
            data.columnInfo.withCheck ~= 0 && ...
            isfield(data, 'checks')
            delete(data.checks(:));
            data.checks = [];
		end
         % restore orig size
		set(data.axes, 'Units', 'pixels');
		pos_ax = get(data.axes, 'position');
		if data.doCheck % restore on left for checkboxes column
            pos_ax(1) = pos_ax(1) - data.checkdx;
            pos_ax(3) = pos_ax(3) + data.checkdx;
		end
        pos_ax(3) = pos_ax(3) + 10; % width of slider;
		set(data.axes, 'position', pos_ax);
		set(data.axes, 'Units', 'normalized');
        delete(data.slider);
        delete(data.btnAdd);
        delete(data.btnDel);
        data = [];
		set(hObj, 'UserData', data);
    end
return

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = editCell(fig, hObj)
% put an edit control over the selected cell
%-------------------------------------------------------------------------
persistent lasttime;

data = get(hObj, 'UserData');
if isempty(data)
    return
end

if isempty(lasttime)
    tic;
    lasttime = toc - .400;
end
thistime = toc;
if (thistime - lasttime) < .350
    mltable(fig, hObj, 'OnDblClick');
    return
end
lasttime = thistime;

set(hObj, 'units', 'pixels');
pt = get(hObj, 'CurrentPoint');
pt=pt(1,:); % strip out 2nd axis info
pos_ax = get(hObj, 'position');
pt(1) = pos_ax(1) + (pt(1) * pos_ax(3));
pt(2) = pt(2) .* pos_ax(4);

d_x = pos_ax(3)/sum(data.columnInfo.weight);
d_y = pos_ax(4)/(data.numRows+1);
  
% find column index
col = -1;
p1 = 0;
p2 = 0;
for i = 1:data.numCols
   p2 = p1 + d_x*data.columnInfo.weight(i);
   if((p1 <= (pt(1)-pos_ax(1))) & (p2 >= (pt(1)-pos_ax(1))))
     col = i;
     break;
   else 
     p1 = p2;
   end;
end
if(col == -1)
  set(hObj, 'units', 'normalized');
  return;
end  

% find row index
row = data.numRows - (floor(pt(2) / d_y));

if(row < 1) % could be header row
  set(hObj, 'units', 'normalized');
  return;
end  

data.selectedCol = col;
data.selectedRow = row + data.ind0;

if isfield(data, 'editBox') && ishandle(data.editBox)
    delete(data.editBox);
end
    
data.hpatch = remakepatch(data.selectedRow, data, hObj);

% continue only if editable    
if(0 ~= data.columnInfo.isEditable(col) && row <= data.numRows)
	unit_d_h = 1/pos_ax(3);
    unit_d_v = 1/pos_ax(4);
	% handle numeric (or not) data
	if data.columnInfo.isNumeric(col)
        ebtxt = data.data{row + data.ind0, col}/data.columnInfo.multipliers(col);
        ebtxt = num2str(ebtxt, data.columnInfo.formats{col});
	else
        ebtxt = UnStripChars(data.data{row + data.ind0, col});
	end
	% set the edt control contents and position
	% callback for entering cell data
	endfcn = sprintf('mltable(%14.13f, %14.13f, ''SetCellValue'');',fig, hObj);
	data.editBox = uicontrol('style', 'edit', 'units', 'pixels',...
        'Callback', endfcn);
	set(data.editBox, 'FontSize', data.gFont.size, 'FontName', data.gFont.name,...
        'FontWeight', 'normal');
	set(data.editBox, 'string', ebtxt, 'UserData', [row col]);
	ext_eb = get(data.editBox, 'extent');
	ext_eb(4) = d_y + 3;
	pos = [(pos_ax(1)-unit_d_h+p1) ,...
            pos_ax(2) + ...                             % start of table 
              ceil((data.numRows - row) * d_y) + ...    % cvt index to row #, get offset from ystart
                 (ceil(d_y) - ext_eb(4))/2, ...          % add half of the ctrl height??
            floor(d_x*data.columnInfo.weight(col))-unit_d_h, ...
            ext_eb(4)];
	% fprintf(1, 'Click at point (%3.2f, %3.2f) on cell (%d, %d) coordinates [%3.2f %3.2f %3.2f %3.2f]\n', ...
	%            pt(1), pt(2), col, row, pos(1), pos(2), pos(3), pos(4));
	set(data.editBox, 'Position', pos, 'HorizontalAlignment' ,'Left');
	set(fig, 'CurrentObject', data.editBox);
	%set(data.editBox,'Editing','true');
end
set(hObj, 'UserData', data);
set(hObj, 'units', 'normalized');
return

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = setCellValue(fig, hObj)
% when edit control calls back, update data in cell
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if isempty(data)
    return
end

if ~isfield(data,'editBox') || ~ishandle(data.editBox)
    return
end

ind = get(data.editBox, 'UserData');
if isempty(ind)
    return;
end;

nums = StripChars(get(data.editBox, 'string'));
row = ind(1) + data.ind0; col = ind(2);
d_old = data.data{row, col};
if data.columnInfo.isNumeric(col)
	num = sscanf(nums, '%f');
	if(isempty(num))
      errordlg('Please enter a valid number', 'Error', 'modal');
      return;
	end     
    d_new = num*data.columnInfo.multipliers(col);
	if(d_old == d_new)
       delete(data.editBox);
       return;
	end
    if mltable(fig, hObj, 'OnSetCell', [], [], [], [], d_old, d_new(1))
        nums = num2str(num, data.columnInfo.formats{col});
        data.data{row, col} = d_new;
    else
        nums = num2str(d_old/data.columnInfo.multipliers(col), data.columnInfo.formats{col});
    end
else
    if mltable(fig, hObj, 'OnSetCell', [], [], [], [], data.data{row, col}, nums)
        data.data{row, col} = nums;
    else
        nums = data.data{row, col};
    end
end

% need to check handles, since OnSetCell callback may have refreshed table
if ishandle(data.editBox)
    delete(data.editBox);
end
if ishandle(data.txtCells(row - data.ind0, col))
    set(data.txtCells(row - data.ind0, col), 'string', nums);
    set(hObj,'UserData',data);
end

return

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = addRow(fig, hObj)
% insert a row into the table
%-------------------------------------------------------------------------
data = get(hObj, 'UserData');
if isempty(data)
    return
end

% increase the number of rows
data.maxRows = data.maxRows + 1;
selRow = data.selectedRow;
if selRow < 1
    selRow = 1;
end

% move data from new row and following down one row
dtmp = data.data;
ctmp = data.isChecked;

for jj = size(dtmp, 1):-1:selRow % for subsequent rows
    ctmp(jj+1) = ctmp(jj);
    for ii = 1:size(dtmp,2) % for all columns
        dtmp{jj+1, ii} = data.data{jj, ii};
    end
end
% zap data in new row
ctmp(selRow) = 0;
for ii = 1:size(dtmp,2) % for all columns
    dtmp{selRow, ii} = '';
end

data.data = dtmp;
data.isChecked = ctmp;

set(hObj, 'UserData', data);
resizeTable(fig, hObj);
return

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function data = delRow(fig, hObj)
% delete the currently selected row
%-------------------------------------------------------------------------
global MINROWS

data = get(hObj, 'UserData');
if isempty(data)
    return
end

selRow = data.selectedRow;
if selRow < 1
    selRow = data.maxRows;
end

if data.maxRows < 2
    return
end

% decrease the number of rows
data.maxRows = data.maxRows - 1;
mr = data.maxRows;
cols = size(data.data,2);
dtmp = {};
ctmp = [];
% remove the selected row
% copy previous data
for ii = 1:selRow-1
    ctmp(ii) = data.isChecked(ii);
	for jj = 1:cols
        dtmp{ii, jj} = data.data{ii, jj};
	end
end    
% copy following data
for ii = selRow+1:size(data.data,1)
    ctmp(ii-1) = data.isChecked(ii);
	for jj = 1:cols
        dtmp{ii-1, jj} = data.data{ii, jj};
	end
end

if data.maxRows < MINROWS
    blanks = cell(1, size(dtmp, 2));
    for ii = data.maxRows+1:MINROWS
        dtmp = [dtmp; blanks];
        ctmp(ii) = 0;
    end
    data.maxRows = MINROWS;
end

data.data = dtmp;
data.isChecked = ctmp;

set(hObj, 'UserData', data);
resizeTable(fig, hObj);
return


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [] = setOnSetCell(fig, hObj, varargin);
% establishes an action for double-clicking on a cell 
% call as "mltable(fig, hobj, 'SetDblClick', 'Myfunc', 'fmt_str', args..."
% where "fmt_str" and args are optional
% function name supplied will be called as
% [] = Myfunc(hObject, [], handles, arg1, arg2, ...)
% ------------------------------------------------------------------------
  sfig  = varargin{1};
  sfunc = varargin{2};
  data = get(hObj, 'UserData');
  data.OnSetCellFcn{1} = sfig;
  data.OnSetCellFcn{2} = sfunc;
%  if nargin < -3
  if nargin > 4
      data.SetCellFmt = varargin{3};
  end
  set(hObj, 'UserData', data);
return

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [res] = onSetCell(fig, hObj, varargin);
% establishes an action for finising a cell edit 
% ------------------------------------------------------------------------
    res = true;
	data = get(hObj, 'UserData');
    if ~isfield(data,'OnSetCellFcn')
        return
    end
    % start with func name, hObject (fig), evdata ([]), and handles
    sfmt  = '%s(''%s'', %14.13f, [], handles';
    sfig  = data.OnSetCellFcn{1};
    sfunc = data.OnSetCellFcn{2};
    if isfield(data,'SetCellFmt')
        sfmt = [sfmt ', ' data.SetCellFmt];
    end
    handles = guidata(fig);
    fcn = sprintf(sfmt, sfig, sfunc, hObj, varargin{:});
    fcn = [fcn ');'];
    try
        res = eval(fcn);
    catch
        res = [];
        disp(['error in OnSetCell callback:\n' lasterr]);
    end
return

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [] = setDblClick(fig, hObj, varargin);
% establishes an action for double-clicking on a cell 
% call as "mltable(fig, hobj, 'SetDblClick', 'Myfunc', 'fmt_str', args..."
% where "fmt_str" and args are optional
% function name supplied will be called as
% [] = Myfunc(hObject, [], handles, arg1, arg2, ...)
% ------------------------------------------------------------------------
  sfig  = varargin{1};
  sfunc = varargin{2};
  data = get(hObj, 'UserData');
  data.OnDblClickFcn{1} = sfig;
  data.OnDblClickFcn{2} = sfunc;
  if nargin < -3
      data.DblClickFmt = varargin{2};      
  end
  set(hObj, 'UserData', data);
return

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [] = onDblClick(fig, hObj, varargin);
% establishes an action for double-clicking on a cell 
% ------------------------------------------------------------------------
	data = get(hObj, 'UserData');
    if ~isfield(data,'OnDblClickFcn')
        return
    end
    % start with func name, hObject (fig), evdata ([]), and handles
    sfmt  = '%s(''%s'', %14.13f, [], handles';
    sfig  = data.OnDblClickFcn{1};
    sfunc = data.OnDblClickFcn{2};
    if isfield(data,'DblClickFmt')
        sfmt = [sfmt ', ' data.DblClickFmt];
    end
    sfmt = [sfmt ');'];
    handles = guidata(fig);
    fcn = sprintf(sfmt, sfig, sfunc, hObj, varargin{:});
    try
        eval(fcn);
    end
return

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [] = onCheck(fig, hObj, row);
% What happens when a user clicks on a row's checkbox
% ------------------------------------------------------------------------
	data = get(hObj, 'UserData');
    rowchk = row + data.ind0;
    vchk = get(data.checks(row),'value');
    data.isChecked(rowchk) = vchk;
    set(hObj,'UserData',data);
return

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [] = setCheck(fig, hObj, row, vchk);
% allows parent to set the data row's checked state
% ------------------------------------------------------------------------
	data = get(hObj, 'UserData');
    if row > data.maxRows
        return
    end
    if row - data.ind0 <= data.numRows
        set(data.checks(row - data.ind0),'value', vchk);
    end
    data.isChecked(row) = vchk;
    set(hObj,'UserData',data);
return


% ------------------------------------------------------------------------
% --------------------------------------------------------------------
function [sout] = StripChars(sin)
    sbad = '_';
    sout = sin;
    stmp = [];
    [sf, sr] = strtok(sin, sbad);
    while ~isempty(sr) % found a '_' char
        stmp = [stmp sf '\_']; % replace with '\_'
        [sf, sr] = strtok(sr, sbad);
    end
    if ~isempty(stmp)
        stmp = [stmp sf];
        sout = stmp;
    end
return

% --------------------------------------------------------------------
function [sout] = UnStripChars(sin)
    sbad = '\\';
    sout = sin;
    stmp = [];
    [sf, sr] = strtok(sin, sbad);
    while ~isempty(sr) % found a '\' char
        stmp = [stmp sf ]; % remove '\'
        [sf, sr] = strtok(sr, sbad);
    end
    if ~isempty(stmp)
        stmp = [stmp sf];
        sout = stmp;
    end
return

function result = trim(string)
[nR, nC] = size(string);

indexStart = 1;
indexEnd   = nC;

for i = 1:nC
  if(string(1, indexStart) == ' ')
    indexStart = indexStart + 1;
  else 
    break;
  end
end  

for i = nC:-1:1
  if(string(1, indexEnd) == ' ')
    indexEnd = indexEnd - 1;
  else 
    break;
  end
end  

result = string(indexStart:indexEnd);

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [data] = mltest
% function to test the table
% ------------------------------------------------------------------------
	fig = figure('renderer','zbuffer');
    pos = get(fig, 'position');
    pos(3) = 440;pos(4)=160;
    set(fig,'position',pos);
    
	tbl = axes('units', 'pixels','position', [20 20 400 120]);
	cell_data = {... 
              'Alpha',   1, 2, 3,'';...
              'Bravo',   4, 5, 6,'';...
              'Charlie', 7, 8, 9,'';...
              'Dog',    10,11,12,'';...
              'Echo',   13,14,15,'';...
              'Foxtrot',16,17,18,'';...
              'Golf',   19,20,21,'';...
              'Hotel',  26,27,28,'';...
              };
	
	columninfo.titles={'Param','Lower Limit','Upper Limit','Initial Value','Result'};
	columninfo.formats = {'%4.6g','%4.6g','%4.6g','%4.6g', '%4.6g'};
	columninfo.weight =      [ 1, .85, .85, .85, .85];
	columninfo.multipliers = [ 1, 1, 1, 1, 1];
	columninfo.isEditable =  [ 1, 1, 1, 1, 0];
	columninfo.isNumeric =   [ 0  1, 1, 1, 1];
	columninfo.withCheck = true;
	columninfo.chkLabel = 'Use';
	rowHeight = 16;
	gFont.size=9;
	gFont.name='Helvetica';
	
	data = mltable(fig, tbl, 'CreateTable', columninfo, rowHeight, cell_data, gFont);
return


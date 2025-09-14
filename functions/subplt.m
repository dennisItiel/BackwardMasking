function theAxis = subplt(nrows, ncols, thisPlot)
% App version : MyPrecious 0.9 (release)
%
%SUBPLT A variation on SUBPLOT with variable spacing between plots.
%	This is a copy of SUBPLOT with the following changes:
%
%	1. The global variable SUBP_FRACT determines the fraction of the total
%	plot size devoted to the axes (default is 0.85). 
%	SUBP_FRACT can be a 2-vector for [horiz vert].
%
%	2. The global variable SUBP_TOTAL determines the fraction of the figure
%	size devoted to the plots (default is 1). 
%	Can also be a 2-vector for [horiz vert].
%
%	3. The global variable SUBP_OFFSET (default 0) can be used to offset
%	the array of plots from their default position (centered on the fig)
%	Can also be a 2-vector for [horiz vert]. Units are fration of fig size.
%
%	4. Special treatments for small fractions or for small number of plots
%	was removed.
%	
%	For more help see SUBPLOT
%
%	Hamutal, April 2005.
%	Last revision: 01 April 2005.
%

global SUBP_FRACT SUBP_TOTAL SUBP_OFFSET
fract = SUBP_FRACT; total = SUBP_TOTAL; offset = SUBP_OFFSET;

if length(fract) == 0, fract = .85; end;
if length(fract) == 1, fract = [fract fract]; end;

if length(total) == 0; total = 1; end;
if length(total) == 1; total = [total total]; end;

if length(offset) == 0; offset = 0; end;
if length(offset) == 1; offset = [offset offset]; end;
offset = offset + (1-total)/2; 

% we will kill all overlapping siblings if we encounter the mnp
% specifier, else we won't bother to check:
narg = nargin;
kill_siblings = 0;
create_axis = 1;
delay_destroy = 0;
if narg == 0 % make compatible with 3.5, i.e. subplot == subplot(111)
	nrows = 111;
	narg = 1;
end

%check for encoded format
handle = '';
position = '';
if narg == 1
	% The argument could be one of 3 things:
	% 1) a 3-digit number 100 < num < 1000, of the format mnp
	% 2) a 3-character string containing a number as above
	% 3) an axis handle
	code = nrows;

	% turn string into a number:
	if(isstr(code)) code = eval(code); end

	% number with a fractional part can only be an identifier:
	if(rem(code,1) > 0)
		handle = code;
		if ~strcmp(get(handle,'type'),'axes')
			error('Requires valid axes handle for input.')
		end
		create_axis = 0;
	% all other numbers will be converted to mnp format:
	else
		thisPlot = rem(code, 10);
		ncols = rem( fix(code-thisPlot)/10,10);
		nrows = fix(code/100);
		if nrows*ncols < thisPlot
			error('Index exceeds number of subplots.');
		end
		kill_siblings = 1;
		delay_destroy = (code == 111);
	end
elseif narg == 2
	% The arguments MUST be the string 'position' and a 4-element vector:
	if(strcmp(lower(nrows), 'position'))
		pos_size = size(ncols);
		if(pos_size(1) * pos_size(2) == 4)
			position = ncols;
		else
			error(['subplot(''position'',',...
			' [left bottom width height]) is what works'])
		end
	else
		error('Unknown command option')
	end
elseif narg == 3
	% passed in subplot(m,n,p) -- we should kill overlaps
	% here too:
	kill_siblings = 1;
end

% if we recovered an identifier earlier, use it:
if(~isempty(handle))
	set(get(0,'CurrentFigure'),'CurrentAxes',handle);
% if we haven't recovered position yet, generate it from mnp info:
elseif(isempty(position))
	if (thisPlot < 1)
		error('Illegal plot number.')
        elseif (thisPlot > ncols*nrows)
                error('Index exceeds number of subplots.')
	else
		row = (nrows-1) -fix((thisPlot-1)/ncols);
		col = rem (thisPlot-1, ncols);

% For this to work the default axes position must be in normalized coordinates
		def_pos = [offset total];
		totalwidth = def_pos(3);
		totalheight = def_pos(4);
		width = (totalwidth * fract(1)) / ncols;	
		height = (totalheight * fract(2)) / nrows;
		position = [def_pos(1)+(col+.5)*totalwidth/ncols-width/2 ...
			    def_pos(2)+(row+.5)*totalheight/nrows-height/2 ...
			    width height];
	end
end

% kill overlapping siblings if mnp specifier was used:
nextstate = get(gcf,'nextplot');
if strcmp(nextstate,'replace'), nextstate = 'add'; end
if(kill_siblings)
	if delay_destroy, set(gcf,'NextPlot','replace'); return, end
	sibs = get(gcf, 'Children');
	for i = 1:length(sibs)
		if(strcmp(get(sibs(i),'Type'),'axes'))
			units = get(sibs(i),'Units');
			set(sibs(i),'Units','normalized')
			sibpos = get(sibs(i),'Position');
			set(sibs(i),'Units',units);
			intersect = 1;
			if(     (position(1) >= sibpos(1) + sibpos(3)) | ...
		                (sibpos(1) >= position(1) + position(3)) | ...
               			(position(2) >= sibpos(2) + sibpos(4)) | ...
		                (sibpos(2) >= position(2) + position(4)))
	               		 intersect = 0;
			end
			if intersect
				if any(sibpos ~= position)
					delete(sibs(i));
				else
					set(gcf,'CurrentAxes',sibs(i));
                                        if strcmp(nextstate,'new')
						create_axis = 1;
					else
						create_axis = 0;
					end
				end
			end
		end
	end
	set(gcf,'NextPlot',nextstate);
end

% create the axis:
if create_axis
	if strcmp(nextstate,'new'), figure, end
	ax = axes('units','normal','Position', position);
        set(ax,'units',get(gcf,'defaultaxesunits'))
else 
	ax = gca; 
end


% return identifier, if requested:
if(nargout > 0)
	theAxis = ax;
end

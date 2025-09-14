function [handle, low, high] = img(matrix, fig_num, low, high)
%
% App version : MyPrecious 0.9 (release)
%
%IMG	Display an image.
%	[HANDLE, LOW, HIGH] = IMG(MATRIX, FIG_NUM, LOW, HIGH)
%	    Hamutal April 2005
%       Last revision: 18 April 2017.

if (nargin == 1); fig_num = gcf; end;
if (nargin == 3); high = low; low = fig_num; fig_num = gcf; end;
if (nargin < 3),
	low  = min(min(matrix));
	high = max(max(matrix));
end;
if (high < low),  disp('warning: high clip is lower than low clip'); end;
if (high == low),
	disp('warning: high clip equals low clip');
	low = low - 1; high = high + 1;		% make it look gray
end;

% Dec 2017: Revised to work with newer (post 2014) Matlab versions. R.O.
if ~isa(fig_num,'double')
    fig_num=fig_num.Number;
end

if fig_num ~= gcf,
	fh = figure(fig_num);
else
	fh = gcf;
end;

if ~isa(fh,'double')
    fh=fh.Number;
end

if ( ~strcmp( get(fh, 'UserData') ,'img_figure') ),
	ori = orient;
	clf reset;
	orient(ori);
	%%colormenu;
%     	uimenu(zoomh, 'label', 'Zoom', 'callback', 'zoom'); % Zoom menu was removed because it wasn't used R.O 2017
% 	uimenu(zoomh, 'label', 'Zoom xy', 'callback', 'zoomxy');
% 	uimenu(zoomh, 'label', 'Zoom x', 'callback', 'zoomx');
% 	uimenu(zoomh, 'label', 'Zoom y', 'callback', 'zoomy');
% 	uimenu(zoomh, 'label', 'Zoom out', 'callback', 'zoomout');
% 	quith = uimenu(fh,'label','Quit','callback',['delete(' num2str(fh) ')']);
% 	set(quith,'BackgroundColor','red','Position',5);
	colormap gray;
	set(fh, 'UserData', 'img_figure');
end;

ncolors = size(colormap, 1);
[m,n] = size(matrix);
matrix = (ncolors-1) * (matrix - low) / (high - low) + 1;

image(matrix);
if ((nargout == 1) | (nargout == 3)); handle = fh; end;
if (nargout == 2); handle = low; low = high; end;

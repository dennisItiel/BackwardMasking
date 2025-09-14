function [handles, low, high] = mimg2(x, xsize, ysize, low, high,    dataStr,  itayColorVec, roi, width)
%MIMG  "Multi" img.
%       [HANDLES, LOW, HIGH] = mimg(X, XSIZE, YSIZE, LOW, HIGH, WIDTH)
%	Takes a set of images in columns and displays them in array of SUBPLTs.
%	It assumes data within each column is ordered by rows
%	(like iv images or experimental data).
%
%	Inputs:
%
%	X	     : the data matrix with an image in each column
%	XSIZE, YSIZE : size of each image. Default: 64 by 64.
%	LOW, HIGH    : Clip values. (default is 'auto').
%		       LOW and HIGH can each be either a scalar or a vector of 
%			   different clip values for each image.
%		       If LOW is 'auto' => autoclip each image separately.
%		       If LOW is 'all'  => use global min/max for all images.
%	WIDTH        : optional width of the result in images. 
%		       Default is ceil(sqrt(# of images)) (same if WIDTH is 0)
%	Can also use MIMG(X, XSIZE, YSIZE, WIDTH). Clipping would be 'auto'.
%
%	Outputs:
%
%	HANDLES	  : Handles for the axes.
%	LOW, HIGH : Vectors of actual clip values used.
%	If only 2 output are asked then LOW, HIGH are supplied.
%	
%	See also: MRESHAPE, IMG, SUBPLT
%
%	Hamutal April 2005.
%	Last revision: 01 April 2005.

if (min(size(x)) <= 1); x = x(:); end;
[m,n] = size(x);
if (nargin < 3); xsize = 64; ysize = 64; end;
if (nargin < 4), low = 'auto'; end;

if (nargin < 9) ; width   = ceil(sqrt(n)); end;
if (width == 0) ; width   = ceil(sqrt(n)); end;
if (nargin == 4) & ~isstr(low); width = low; low = 'auto'; end;
if ~exist('dataStr','var'); dataStr = zeros(1,n); end;
if ~exist('itayColorVec','var'); itayColorVec = [0 0 0]; end;
if ~exist('roi','var');
    b = []; 
else
    b = zeros(10000,1); b(roi) = 1;
end;

if itayColorVec == 'b'
    itayColorVec = [0 0 0];
elseif itayColorVec == 'w'
    itayColorVec = [1 1 1];
end


xsize = round(xsize);
ysize = round(ysize);

auto = 0;
if isstr(low),
	low = lower(low);
	if strcmp(low, 'auto'),		% auto clip each image separately
		auto = 1;
	elseif strcmp(low, 'all'),	% global min and max for all
		low  = min(min(x));
		high = max(max(x));
	else
		error('String argument LOW should be ''all'' or ''auto''!');
	end;
end;
if ~auto,
	if length(low) == 1,
		low = low(ones(1,n));
	end;
	if length(high) == 1,
		high = high(ones(1,n));
	end;
	if (length(low) ~= n) | (length(high) ~= n),
		error('LOW and HIGH should have 1 or n-images elements!')
	end;
end;
height = ceil(n/width);
img(1:10);
subplt(height, width, 1);
axs = []; lows = []; highs = [];
for i=1:n,
	ax = subplt(height, width, i);
	if auto,
	   [l, h] = img(mreshape(x(:,i), 0, 0, xsize, ysize));hold on;
       text(6,10,num2str(dataStr(i)),'color',itayColorVec,'FontWeight','bold')%
       if ~isempty(b)
           contour(reshape(b,[100,100])','w');
       end
    else
	   [l, h] = img(mreshape(x(:,i), 0, 0, xsize, ysize), low(i), high(i));hold on;
       text(6,10,num2str(dataStr(i)),'color',itayColorVec,'FontWeight','bold')%
       if ~isempty(b)
           contour(reshape(b,[100,100])','w');
       end;
    end;
	axis image; axis off;
	axs   = [axs; ax];
	lows  = [lows; l];
	highs = [highs; h];
end;

low = lows;
high = highs;

if ((nargout == 1) | (nargout == 3)); handles = axs; end;
if (nargout == 2); handles = low; low = high; end;

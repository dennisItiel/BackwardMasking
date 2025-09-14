function y = mreshape(x, width, border, xsize, ysize)
%
% App version : MyPrecious 0.9 (release)
%
%MRESHAPE  "Multi" reshape.
%       Y = mreshape(X, WIDTH, BORDER, XSIZE, YSIZE)
%	Takes a set of images in columns and arranges them in an array.
%	It assumes data within each column is ordered by rows
%	(like iv images or experimental data).
%	X: the data matrix with an image in each column
%	WIDTH: optional width of the result in images. 
%	       Default: ceil(sqrt(# of images)) (sets to this also if eq 0)
%	BORDER: optional # of pixels inserted between images and filled with 
%	        min(X). Default: 0 if x is a vector, 1 otherwise.
%		Can also have a second element which determines the value at
%		the border. [0,1] corresponds to [min(X),max(X)]
%	XSIZE, YSIZE: size of each image. Default: 64 by 64.
%	
%	Similar (but simpler) to NEW2D.
%
%	Hamutal Slovin, April 2005.
%	Last revision: 01 April 2005.
%

if (min(size(x)) <= 1); x = x(:); end;
[m,n] = size(x);
if (nargin < 5); xsize = 64; ysize = 64; end;
if (nargin < 3); border  = (n>1); end;
if (nargin < 2) ; width   = ceil(sqrt(n)); end;
if (width == 0) ; width   = ceil(sqrt(n)); end;

minx = min(min(x));
if length(border) == 2,
	maxx = max(max(x));
	back_color = minx + border(2) * (maxx - minx);
else,
	back_color = minx; 
end;
border = border(1);

xsize = floor(xsize);
ysize = floor(ysize);

xstep=xsize+border;
ystep=ysize+border;
y=zeros(ystep*ceil(n/width) + border, xstep*width + border) + back_color;
for i=0:(n-1),
  y([1:ysize] + ystep*floor(i/width) + border, [1:xsize] + xstep*rem(i,width) + border) = reshape(x(:,i+1), xsize, ysize)';
end;

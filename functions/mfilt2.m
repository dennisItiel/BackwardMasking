function y = mfilt2(x, xsize, ysize, s, type, s1)

%MFILT2	"Multi" 2D Gaussian filter
%	Y = MFILT2(X, XSIZE, YSIZE, S, TYPE, S1) 
%
%	X:	Matrix to be filtered with 2D images at each column
%	XSIZE, YSIZE size of each image 
%		This order assumes the images are ordered by rows,
%		otherwise the order should be YSIZE, XSIZE.
%	S:	Sigma of gaussian to use (default 1)
%	TYPE:	A string of 1-2 characters to determind type of filter.
%		'l' (default) for low-pass, 'h' for high-pass,
%		'b' for band-pass (difference of gaussians)
%		If also 'm' is given (e.g, 'lm'), the image will be padded
%		with it's mirror reflections before filtering, instead of the
%		default which effectively does wrap around.
%	S1:	Sigma for high pass in case of band-pass (default 1.5*S)
%
%	See also: FILT2, FILTX, FILTY, MFILTX, MFILTY
%
%	HS April 2005.
%	Last revision: 01 Apr 2005.

if (nargin < 4) 
    s = 1; 
end;
if (nargin < 5) 
    type = 'l'; 
end;
if (nargin < 6) 
    s1 = 1.5 * s; 
end;

if (~ischar(type) || length(type(:)) > 2),
	error('TYPE should be a 1- or 2-element string vector!');
end;

xsize = floor(xsize);
ysize = floor(ysize);
[n, m] = size(x);

if (xsize*ysize ~= n), error('column length is not xsize*ysize'); end;

mind = find(type == 'm');
do_mirror = length(mind);
if do_mirror,
	type(mind) = [];
	pad = ceil(4*s);	% Padding with 4*sigma pixels.
	if (type == 'b'), pad = ceil(4*max(s,s1)); end;
	pm = min(pad, ysize);	% Pad is not more than image size
	pn = min(pad, xsize);
	ysize_orig = ysize;
	xsize_orig = xsize;
	ysize = ysize + 2*pm;
	xsize = xsize + 2*pn;
end;

sf  = 1 / (2*pi*s);	% sigma in frequency space
sf1 = 1 / (2*pi*s1);	

% xsize and ysize are in "wrong" order so I can use reshape and not mreshape
%[xx,yy]= meshgrid(fftshift([-ysize/2:1:(ysize/2-1)]/ysize), fftshift([-xsize/2:1:(xsize/2-1)]/xsize));
% The above is correct only for even m and n. The next 3 lines are more general

yind = [0:ceil(ysize/2-1), ceil(-ysize/2):-1] / ysize;
xind = [0:ceil(xsize/2-1), ceil(-xsize/2):-1] / xsize;
[yy,xx]= meshgrid(yind, xind);

g = exp(-(xx.^2+yy.^2)/(2*sf.^2));
if (type =='h'), g = 1-g; end;
if (type =='b'),
	g1 = exp(-(xx.^2+yy.^2)/(2*sf1.^2));
	g = g - g1;
end;

y = zeros(size(x));
if ~do_mirror,
	for i=1:m,
		fx = fft2(reshape(x(:,i), xsize, ysize));
		tmp = real(ifft2(fx.*g));
		y(:,i) = tmp(:);
	end;	
else
		ypad_ind = [pm:-1:1 1:ysize_orig ysize_orig:-1:(ysize_orig-pm+1)];
		xpad_ind = [pn:-1:1 1:xsize_orig xsize_orig:-1:(xsize_orig-pn+1)];
	for i=1:m,
		fx = reshape(x(:,i), xsize_orig, ysize_orig);
		fx = fx(xpad_ind, ypad_ind);
		fx = fft2(fx);
		tmp = real(ifft2(fx.*g));
		tmp = tmp([1:xsize_orig]+pn, [1:ysize_orig]+pm);	% get back to original size
		y(:,i) = tmp(:);
	end;	
end;

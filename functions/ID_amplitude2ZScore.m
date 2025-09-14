function [VSD_zs] = ID_amplitude2ZScore(VSD_matfile)
%UNTITLE 
%   Convert amplitude to Z-Score from baseline

BASE_LINE = 10:23;

if ndims(VSD_matfile) == 3
    [~,mu,sigma] = zscore(VSD_matfile(:,BASE_LINE,:), 0, 2);
elseif ndims(VSD_matfile) == 2
    [~,mu,sigma] = zscore(VSD_matfile(:,BASE_LINE), 0, 2);
end

VSD_zs =(VSD_matfile-mu)./sigma;

end
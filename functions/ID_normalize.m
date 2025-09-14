function [vertical_data,horizontal_data] = ID_normalize(condsn2_bl,condsn4_bl,frames,v1_pix,normType)
%NORMALIZE before training model
%   INPUT:
%   condsn2_bl = vertical data  10000 X 256 X n
%   condsn4_bl = horizontal data    10000 X 256 X n
%   normType = type of normalization
%
%   OUTPUT:
%   vertical_zs, horizontal_zs: n X 256 X 10000

if strcmp(normType, 'zScore')
    normalization = @ ID_amplitude2ZScore;
elseif strcmp(normType, 'iqr')
    normalization = @ ID_amplitude2IQR;
end

% Vertical frames upload
try
    vertical_data = normalization(condsn2_bl(v1_pix,:,:));  % convert to Z score
    vertical_data = permute(vertical_data,[3 2 1]); %change dimensions order
    vertical_data = vertical_data(:,frames,:);
catch
    vertical_data = [];
end

% Horizontal frames upload
try
    horizontal_data = normalization(condsn4_bl(v1_pix,:,:));
    horizontal_data = permute(horizontal_data,[3 2 1]); %change dimensions order
    horizontal_data = horizontal_data(:,frames,:);
catch
    horizontal_data = [];
end
end
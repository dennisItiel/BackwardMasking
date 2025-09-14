function [vertical_in, vertical_out, horizontal_in, horizontal_out] = ID_leaveOneOut(vertical_data, horizontal_data)
%LEAVE ONE OUT
%   separetes data randomly
%   vertical and horizontal 'in' has same size
%   data => trials X frames X 100000

num_ver = size(vertical_data, 1);
num_hor = size(horizontal_data, 1);
num_of_in = min(num_ver, num_hor) - 1; % number of trials for trainig set

% assign vertical data
ver_in_inx = false(num_ver, 1);
ver_in_inx(randperm(num_ver, num_of_in)) = 1;
vertical_in = vertical_data(ver_in_inx,:,:);
vertical_out = vertical_data(~ver_in_inx,:,:);

% assign horizontal data
hor_in_inx = false(num_hor, 1);
hor_in_inx(randperm(num_hor, num_of_in)) = 1;
horizontal_in = horizontal_data(hor_in_inx,:,:);
horizontal_out = horizontal_data(~hor_in_inx,:,:);

end
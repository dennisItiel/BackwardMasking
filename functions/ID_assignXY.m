function [X,Y] = ID_assignXY(vertical_data,horizontal_data)
%ASSIGN data o VSDI to readable format for SVM
%   prepare the X & Y
%   X -> data yo learn from (samples X pixels)
%   Y -> labels vector of 1,-1 (sample X 1)


vertical_data = permute(vertical_data, [2 1 3]);
A = reshape(vertical_data, size(vertical_data,1) * size(vertical_data,2), []);

horizontal_data = permute(horizontal_data, [2 1 3]);
B = reshape(horizontal_data, size(horizontal_data,1) * size(horizontal_data,2), []);

X = [A;B];
Y = [ones(0.5 * size(X,1),1); -ones(0.5 * size(X,1),1)];


% end
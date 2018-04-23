function [X, Y] = arxinput(Z, N)
%ARXINPUT returns the linear regressor input output matrices needed
%   Detailed explanation goes here
na = N(1); % number of past outputs in regressor
nb = N(2); % number of past inputs
nk = N(3); % delays in input signal

[n, d] = size(Z);

% i_start = max(1+na, nk+nb);
% if n - i_start > 0
%     X = zeros(n-i_start, na+nb); % regressor matrix
%     Y = zeros(n-i_start, 1);
%     for i = 1:(n-i_start)
%         Iy = (i_start+i-1):-1:(i_start+i-na);
%         Iu = (i_start+i-1-nk):-1:(i_start+i-nk-nb);
%         X(i, :) = [-Z(Iy,1)', Z(Iu,2)'];
%         Y(i,1) = Z(i+i_start, 1);
%     end
% end


% new with padding

XU = zeros(n+nb-1+nk, nb);
for i = 1:nb
    XU(i+nk+1:i+nk+n, i) = Z(:, 2);
end
XU = XU(1:n,:);

XY = zeros(n+na, na);
for i = 1:na
    XY(i+1:i+n, i) = Z(:, 1);
end
XY = XY(1:n,:);

X = [-XY, XU];
Y = Z(:, 1);
end


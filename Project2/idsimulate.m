function y_sim = idsimulate(m, Z)
%IDSIMULATE Summary of this function goes here
%   Detailed explanation goes here

method = 'build_in';

switch method
    case 'build_in'
        H = id2tf(m);
        y_sim = lsim(H, Z(:,2), (1:size(Z, 1))');
    otherwise
        n = size(Z, 1);
        y_sim = zeros(n, 1);

        u = Z(:,2);
        for i = 1:n
            ui = u(1:i, 1);
            yi = y_sim(1:i, 1);
            zi = [yi, ui];
            [X, ~] = arxinput(zi, m.n);
            yi_sim = X(end,:)*m.theta;
            y_sim(i, 1) = yi_sim;
        end
end
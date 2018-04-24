function y_pred = idpredict(m, z, horizon)

n = size(z, 1);
Y_pred = zeros(n, horizon);
y = z(:, 1); u = z(:,2);

U = zeros(n, m.n(2)); Y = zeros(n, m.n(1));

for i = 1:m.n(2)
    U(i+1+m.n(3):end, i) = u(1:end-i-m.n(3));
end

for i = 1:m.n(1)
    Y(i+1:end, i) = -y(1:end-i);
end

switch m.model
    case 'OE'
        y_pred = idsimulate(m, z);  % prediction of OE is the same as simulation, regardless of horizon!
    case 'ARX'
        for i = 1:horizon
            X = [Y, U];
            Y_pred(:, i) = X*m.theta;
            % update the regressor matrix
            Y(2:end, 2:end) = Y(1:end-1, 1:end-1);
            Y(2:end, 1) = -Y_pred(1:end-1, i);
%             U(2:end, :) = U(1:end-1, :);
%             U(1, :) = zeros(size(U(1, :)));
        end
        y_pred = Y_pred(:, end);
    otherwise
end
  
end
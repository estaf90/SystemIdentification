function y_pred = idpredict(m, z, horizon)

n = size(z, 1);
Y_pred = zeros(n, horizon);
y = z(:, 1); u = z(:,2);

% X = zeros(n, m.n(1)+ m.n(2));
U = zeros(n, m.n(2)); Y = zeros(n, m.n(1));

for i = 1:m.n(2)
    U(i+1+m.n(3):end, i) = u(1:end-i-m.n(3));
end

for i = 1:m.n(1)
    Y(i+1:end, i) = -y(1:end-i);
end

% switch m.model
%     case 'OE'
%         Y_pred(:, i) = X*m.theta;
%     case ''
for i = 1:horizon
    X = [Y, U];
    switch m.model
        case 'ARX'
            Y_pred(:, i) = X*m.theta;
        case 'OE'
            Y_pred(:, i) = X*m.theta;  % wrong...
        case 'OE...'
            Y_pred(:, i) = X(:, m.n(1)+1:end)*m.theta(m.n(1)+1:end);
            if i == 1
                e = y-Y_pred(:, i);
            else
                e = Y_pred(:, i-1)-Y_pred(:, i);
            end
            [E,~] = arxinput([e, e], m.n);
            Y_pred(:, i) = Y_pred(:, i) + E(:, m.n+1:end)*m.theta(m.n(1)+1:end);
        otherwise
            fprintf('idpredict is not implemented for model %s', m.model)
    end
    % update the regressor matrix
    Y(2:end, 2:end) = Y(1:end-1, 1:end-1);
    Y(2:end, 1) = -Y_pred(1:end-1, i);
end

y_pred = Y_pred(:, end);
  
end
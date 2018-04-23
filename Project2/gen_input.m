function u = gen_input(n, mode)
%GEN_INPUT Summary of this function goes here
%   Detailed explanation goes here
N = n;
switch mode
    case 'white'
        u = (rand(N,1)-0.5);  % noise process driving the system
    case 'randi'
        u = randi(2, N, 1)-1;
    case 'steps'
        u = zeros(N,1);
        y = 0;
        for i =1:n
            if mod(i, 50) == 10
                y = randi(5);
            end
            u(i) = y;
        end
end

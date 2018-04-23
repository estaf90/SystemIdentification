function y_gen = gen_oe_data(u, a, b, e)
%GEN_DATA Summary of this function goes here
%   Detailed explanation goes here

na = length(a);
nb = length(b);
u_gen = zeros(1, nb);
y_gen = zeros(1, na);

for i = 1:length(u)
    nu = length(u_gen); ny = length(y_gen);
    idxy = ny:-1:(ny-na+1);
    idxu = nu:-1:(nu-nb+1);
    y_gen = [y_gen, -a*y_gen(idxy)' + b*u_gen(idxu)'];
    u_gen = [u_gen, u(i)];
end

y_gen = y_gen((na+1):end);
y_gen = y_gen + e*randn(size(y_gen));

y_gen = y_gen';
end


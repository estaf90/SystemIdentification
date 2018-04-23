function transfer_function = id2tf(m)
  % Return a MATLAB transfer function tf from my arx model
na = m.n(1);
nb = m.n(2);
nk = m.n(3);
diff = nb+nk-na;  % if positive append denominator, else ...
denominator = [1; m.theta(1:na)];
numerator = m.theta((na+1):end);

if diff < 0
    numerator = [numerator; zeros(diff, 1)];
end
if diff > 0
    denominator = [denominator; zeros(diff, 1)];
end
    
Denominator = {denominator'};
Numerator = {numerator'};


transfer_function = tf(Numerator, Denominator, -1);
end
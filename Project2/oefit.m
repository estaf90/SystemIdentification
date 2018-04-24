function m_oe = oefit(Z, N)
%OEFIT Summary of this function goes here
%   Detailed explanation goes here
nb = N(1); nf = N(2); nk = N(3);
y = Z(:,1); u = Z(:, 2);

m_arx = arxfit(Z, [4*nf, 4*nb, nk]);  % Fit higher order ARX

% simulate the ARX model giving output ys
ys = idsimulate(m_arx, Z);

method = 'approximate';

switch method
    case 'approximate' % then approximate method
        m_oe = arxfit([ys, u], [nf, nb, nk]);
    case 'optimal'
        [X_ys, ~] = arxinput([ys, ys], m_arx.n);
        [X_u, ~] = arxinput([u, y], m_arx.n);
        X_ys = -X_ys; X_u = -X_u;  % due to how arxinput works...

        X_ys = X_ys(:, 1:m_arx.n(1)); X_u = X_u(:, 1:m_arx.n(1));
        ys_fir = ys + X_ys*m_arx.A;
        u_fir = u + X_u*m_arx.A;

        Z_fir = [ys_fir, u_fir];

        m_oe = arxfit(Z_fir, N);
    otherwise
        disp('Method not implemented')
end

m_oe.model = 'OE';
m_oe.label = sprintf('OE [nb=%d, nf=%d, nk=%d]',nb,nf,nk);
end


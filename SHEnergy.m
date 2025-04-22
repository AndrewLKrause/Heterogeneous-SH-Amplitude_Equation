function E = SHEnergy(U, F, r, params)
    r = r(1);
    eps = params.eps; dx = params.dx;
    Lap = @(u)(eps/dx)^2*[u(2)-u(1), u(1:end-2)-2*u(2:end-1)+u(3:end),u(end-1)-u(end)];
    E = trapz(dx*(-(r/2)*(U).^2+(1/2)*(U+Lap(U)).^2-F(U)));
end
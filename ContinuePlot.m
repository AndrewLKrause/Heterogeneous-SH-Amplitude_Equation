function [rc, mp, rvec, Evec, Avec] = ContinuePlot(f, F, params)

N = params.N; eps = params.eps; dx = params.dx; T = params.T;

Lap = @(u)(eps/dx)^2*[u(2)-u(1); u(1:end-2)-2*u(2:end-1)+u(3:end);u(end-1)-u(end)];
Bih = @(u)Lap(Lap(u));
e = ones(N,1); % Vector of ones to use across the diagonals
JP = spdiags([e e e e e], -2:2, N, N); % five-diagonal sparsity pattern

options = optimoptions('fsolve','Algorithm','trust-region','JacobPattern',JP,'display','off');

%x0 is needed to make vectors the right size.
x0 = 0*linspace(0,1,params.N)';
r = 0.1 +x0; qc = 1+x0;
FRHS = @(U)r.*U-((qc.^2).*U+qc.*Lap(U)+Lap(qc.*U)+Bih(U)) +f(U);
U = RunSimulation(r, qc, f, params,1e-2*randn(N,1));U = U(end,:);
r = r-0.09;
U = RunSimulation(r, qc, f, params,U);U = U(end,:);

Us = fsolve(FRHS,U',options);
error = sqrt(trapz(dx*(Us'-U).^2))
U = RunSimulation(r, qc, f, params,Us');U = U(end,:);

E = SHEnergy(U,F,r,params); A = max(U)-min(U);

ETol = 1e-7;

i = 2; Evec = E; rvec = r(1); Avec = A;

while(A > 1e-6)
    r = r-0.001;
    Uold = U;
    U = RunSimulation(r, qc, f, params,U);
    U = U(end,:);
    FRHS = @(U)r.*U-((qc.^2).*U+qc.*Lap(U)+Lap(qc.*U)+Bih(U)) +f(U);
    Us = fsolve(FRHS,U',options);
error = sqrt(trapz(dx*(Us'-U).^2))
    U = RunSimulation(r, qc, f, params,Us');U = U(end,:);
    E = SHEnergy(U,F,r,params); A = max(U)-min(U);

    Evec(i) = E; rvec(i) = r(1); Avec(i) = A;
    i = i+1;
end
Evec(end) = [];
rc = (rvec(end)+rvec(end-1))/2;
mp = rvec(Evec==min(abs(Evec)));

end
function [U] = RunSimulation(r, qc, f, params, Uinit)

N = params.N; eps = params.eps; dx = params.dx; T = params.T;

% %Form the Laplacian
% e = ones(N,1); % Vector of ones to use across the diagonals
% Lap= spdiags([e -2*e e], -1:1, N, N); % Diagonal Laplacian
% %Lap(1,1) = -1; Lap(N,N) = -1; % Neumann boundary conditions
% Lap(N,1) = 1; Lap(1,N) = 1; % Periodic BCs
% Lap = (eps/dx)^2*Lap; % Scale the finite-difference operator

Lap = @(u)(eps/dx)^2*[u(2)-u(1); u(1:end-2)-2*u(2:end-1)+u(3:end);u(end-1)-u(end)];
Bih = @(u)Lap(Lap(u));
e = ones(N,1); % Vector of ones to use across the diagonals
JP = spdiags([e e e e e], -2:2, N, N); % five-diagonal sparsity pattern

F = @(t,U)r.*U-((qc.^2).*U+qc.*Lap(U)+Lap(qc.*U)+Bih(U)) +f(U);

%opts = odeset('JPattern',JP);%,'reltol',1e-9,'AbsTol',1e-9);
opts = odeset('JPattern',JP,'reltol',params.tol,'AbsTol',params.tol);

[~,U] = ode15s(F, T, Uinit,opts);
end
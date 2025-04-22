function [A] = RunGLESimulation(rp, tp, alpha, beta, params, Uinit)

N = params.N; eps = params.eps; dx = params.dx; T = params.T;

x = linspace(0,1,N)';
gamma = alpha/eps^(1/3);

% %Form the Laplacian
% e = ones(N,1); % Vector of ones to use across the diagonals
% Lap= spdiags([e -2*e e], -1:1, N, N); % Diagonal Laplacian
% %Lap(1,1) = -1; Lap(N,N) = -1; % Neumann boundary conditions
% Lap(N,1) = 1; Lap(1,N) = 1; % Periodic BCs
% Lap = (eps/dx)^2*Lap; % Scale the finite-difference operator

D = (4*eps^(4/3))/dx^2;

Lap = @(u)D*[u(2)-u(1); u(1:end-2)-2*u(2:end-1)+u(3:end);u(end-1)-u(end)];

e = ones(N,1); % Vector of ones to use across the diagonals
JP = spdiags([e e e e e], -1:1, N, N); % three-diagonal sparsity pattern

F = @(t,A)Lap(A)+(rp*(x-tp).*A)/(eps^(2/3))+(3/4)*gamma*A.^3-(10/16)*beta*A.^5;

%opts = odeset('JPattern',JP);%,'reltol',1e-9,'AbsTol',1e-9);
opts = odeset('JPattern',JP,'reltol',params.tol,'AbsTol',params.tol);

[~,A] = ode15s(F, T, Uinit,opts);

A = A(end,:)*eps^(1/6);
end
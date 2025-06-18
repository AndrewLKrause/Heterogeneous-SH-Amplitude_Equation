function [R] = RunGLESimulation(rp, x_0, b,c,d,e, params, Uinit)

N = params.N; eps = params.eps; dx = params.dx; T = params.T;

x = linspace(0,1,N)';
%gamma = c/eps^(1/3);
c_0=-38*b^2/27
c_2 = (c-c_0)/eps^(1/3)

% %Form the Laplacian
% e = ones(N,1); % Vector of ones to use across the diagonals
% Lap= spdiags([e -2*e e], -1:1, N, N); % Diagonal Laplacian
% %Lap(1,1) = -1; Lap(N,N) = -1; % Neumann boundary conditions
% Lap(N,1) = 1; Lap(1,N) = 1; % Periodic BCs
% Lap = (eps/dx)^2*Lap; % Scale the finite-difference operator

D = (4*eps^(4/3))/dx^2;

Lap = @(u)D*[u(2)-u(1); u(1:end-2)-2*u(2:end-1)+u(3:end);u(end-1)-u(end)];

%Adv = (1/(2*dx))*spdiags([ones(N,1),-ones(N,1)],[1,-1],N,N);
Adv = @(u)(eps^(2/3)/(2*dx))*[-u(2); -u(1:end-2)+u(3:end);-u(end-1)];

% Neumann boundary conditions
% Lap(1,1) = -1; Lap(end,end) = -1;
% Adv(1,2)=0; Adv(end,end-1)=0;

e_v = ones(N,1); % Vector of ones to use across the diagonals
JP = spdiags([e_v e_v e_v e_v e_v], -1:1, N, N); % three-diagonal sparsity pattern
%JP = sparse([JP, JP; JP JP]);
%RI = (1:N)'; phiI = (N+1:2*N)';

%A = Uinit;
%size(A(RI))
%size(A(phiI))

%size(Lap(A(RI)) - 4*A(RI).*(Adv(A(phiI)).^2) - (32/27)*b^2*A(RI).^3.*Adv(A(phiI)) + (rp*(x-x_0).*A(RI))/(eps^(2/3))...
%    + 3*c_2*A(RI).^3 + (-3920*b^4/81+116*b*d/3+10*e)*A(RI).^5)
%size(Lap(A(phiI)) + 8*Adv(A(RI)).*Adv(A(phiI))./A(RI) + (32/27)*b^2*A(RI).*Adv(A(RI)))
%size([Lap(A(RI)) - 4*A(RI).*(Adv(A(phiI)).^2) - (32/27)*b^2*A(RI).^3.*Adv(A(phiI)) + (rp*(x-x_0).*A(RI))/(eps^(2/3))...
%    + 3*c_2*A(RI).^3 + (-3920*b^4/81+116*b*d/3+10*e)*A(RI).^5;...
%    Lap(A(phiI)) + 8*Adv(A(RI)).*Adv(A(phiI))./A(RI) + (32/27)*b^2*A(RI).*Adv(A(RI))])


%F = @(t,A)Lap(A) + (rp*(x-x_0).*A)/(eps^(2/3)) + 3*c_2*A.^3 + (-3920*b^4/81+116*b*d/3+10*e)*A.^5;
%F = @(t,A)[Lap(A(RI)) - 4*A(RI).*(Adv(A(phiI)).^2) - (32/27)*b^2*A(RI).^3.*Adv(A(phiI)) + (rp*(x-x_0).*A(RI))/(eps^(2/3))...
%    + 3*c_2*A(RI).^3 + (-3920*b^4/81+116*b*d/3+10*e)*A(RI).^5;...
%    (Lap(A(phiI)) + 8*Adv(A(RI)).*Adv(A(phiI))./A(RI) + (32/27)*b^2*A(RI).*Adv(A(RI)))];
%F = @(t,A)[Lap(A) - 4*A.*(Adv(A(phiI)).^2) - (32/27)*b^2*A(RI).^3.*Adv(A(phiI)) + (rp*(x-x_0).*A(RI))/(eps^(2/3))...
%    + 3*c_2*A(RI).^3 + (-3920*b^4/81+116*b*d/3+10*e)*A(RI).^5;

F = @(t,A)Lap(A) + (rp*(x-x_0).*A)/(eps^(2/3)) + 3*c_2*A.^3 + ((32*b^2-4)/(27^2)-3920*b^4/81+116*b*d/3+10*e)*A.^5;


%size(Lap(Uinit))
%size(rp)
%size(x-x_0.*Uinit)
%size(c_2*Uinit.^3)
%size(10*e*Uinit.^5)

%F(0,Uinit)

%opts = odeset('JPattern',JP);%,'reltol',1e-9,'AbsTol',1e-9);
opts = odeset('JPattern',JP,'reltol',params.tol,'AbsTol',params.tol);

[~,A] = ode15s(F, T, Uinit,opts);

%R = 2*A(:,RI)*eps^(1/6);
R = 2*A*eps^(1/6);
%phi = A(:,phiI);
end
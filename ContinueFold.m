function [rc, maxwellpt] = ContinueFold(f, F, params)
TEnd = 1e4;

%x0 is needed to make vectors the right size.
x0 = 0*linspace(0,1,params.N)';
r = 0.1 +x0; qc = 1+x0;
N = params.N; eps = params.eps; dx = params.dx;

Lap = @(u)(eps/dx)^2*[u(2)-u(1); u(1:end-2)-2*u(2:end-1)+u(3:end);u(end-1)-u(end)];
Bih = @(u)Lap(Lap(u));
e = ones(N,1); % Vector of ones to use across the diagonals
JP = spdiags([e e e e e], -2:2, N, N); % five-diagonal sparsity pattern

options = optimoptions('fsolve','Algorithm','trust-region','JacobPattern',JP,'display','off');
FRHS = @(U,r)r.*U-((qc.^2).*U+qc.*Lap(U)+Lap(qc.*U)+Bih(U)) +f(U);



params.T = [0,TEnd];
U = RunSimulation(r, qc, f, params,1e-1*randn(params.N,1));U = U(end,:);
Us = fsolve(@(U)FRHS(U,r),U',options);
U = RunSimulation(r, qc, f, params,Us');U = U(end,:);
r = r-0.05;
U = RunSimulation(r, qc, f, params,U);U = U(end,:);
Us = fsolve(@(U)FRHS(U,r),U',options);
U = RunSimulation(r, qc, f, params,Us');U = U(end,:);

params.T = [0,TEnd];
ETol = 1e-7;
maxwellpt=0;
while(max(U)-min(U) > 1e-5)
    r = r-0.05; r(1)
    Uold = U;
    U = RunSimulation(r, qc, f, params,U);U = U(end,:);
    Us = fsolve(@(U)FRHS(U,r),U',options);
    U = RunSimulation(r, qc, f, params,Us');U = U(end,:);
    E = SHEnergy(U,F,r,params)
    if(E>ETol && maxwellpt==0)
        re = r+0.05; re(1)
        Ue = Uold;
        E = SHEnergy(Ue,F,re,params)
        while(E<ETol)
            re = re-0.01; re(1)
            Ueold = Ue;
            Ue = RunSimulation(re, qc, f, params,Ue);Ue = Ue(end,:);
            Us = fsolve(@(U)FRHS(U,re),Ue',options);
            Ue = RunSimulation(re, qc, f, params,Us');Ue = Ue(end,:);
            E = SHEnergy(Ue,F,re,params)
        end
        re = re+0.01; Ue = Ueold;
        E = SHEnergy(Ue,F,re,params);
        while(E<ETol)
            re = re-0.001; re(1)
            Ue = RunSimulation(re, qc, f, params,Ue);Ue = Ue(end,:);
            Us = fsolve(@(U)FRHS(U,re),Ue',options);
            Ue = RunSimulation(re, qc, f, params,Us');Ue = Ue(end,:);
            E = SHEnergy(Ue,F,re,params)
        end
        maxwellpt = re(1)+0.0005;
    end
end
r = r+0.05; U = Uold;
while(max(U)-min(U) > 1e-5)
    r = r-0.01; r(1);
    Uold = U;
    U = RunSimulation(r, qc, f, params,U);U = U(end,:);
    Us = fsolve(@(U)FRHS(U,r),U',options);
    U = RunSimulation(r, qc, f, params,Us');U = U(end,:);
end
r = r+0.01; U = Uold;
while(max(U)-min(U) > 1e-5)
    r = r-0.001; r(1);
    U = RunSimulation(r, qc, f, params,U);U = U(end,:);
    Us = fsolve(@(U)FRHS(U,r),U',options);
    U = RunSimulation(r, qc, f, params,Us');U = U(end,:);
end
rc = r(1)+0.0005;
end

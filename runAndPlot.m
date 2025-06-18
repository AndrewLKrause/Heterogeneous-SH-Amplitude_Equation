%This script runs and plots simulations in 1D. Depending on the flags it
%will call continuation code to produce the lines in the figures
%corresponding to local folds in a subcriticfal instability, as well as the
%Maxwell point. It now also computes a steady state amplitude and plts it
%alongise the solution.

clear;rng(3);

%Note that epsilon will be squared in the Laplacian (and taken to the fourth power in the biharmonic) 
%so can become very small numerically.
eps = 2e-3;


b = -2; c = -5.6; d = 1; e = -2;

f = @(U)b*U.^2+c*U.^3+d*U.^4+e*U.^5;

%Capital F is the formal integral of little f, and is only needed if
%subcritical=1.
%F = @(U)(alpha/4)*U.^4-(beta/6)*U.^6;

%Number of grid points and time interval to simulate over.
N = 2000; T = [linspace(0,50000,1e3)];
dx = 1/(N-1); % Spacing between grid points
params.dx = dx; params.N = N; params.eps = eps; params.T = T;
params.tol = 1e-6; 
x = linspace(0,1,N)'; 

% rf is used as an analytical representation of r(x), whereas r is a
% numerical vector evaluated on points in [0,1].
rf = @(x)-cos(pi*x);%+0.3*sin(t*eps^(4/3));
r = rf(x);

%For the amplitude equation, we need the derivative of r evaluated at the
%Turing point. This could be done more automatically if desired but for now
%it needs to be provided.
x_0 = fzero(rf,0.5);
%sp = @(t)0.3*eps^(4/3)*cos(t*eps^(4/3));
rp = pi;


%Either load the initial condition from a file, or use a small
%normally-distributed random one.
%load('init3.mat');load('init_eps_2e-3.mat');
Us = 1e-1*randn(params.N,1);
%load('init_eps_2e-4.mat');


% Run the simulation.
U = RunSimulation(r, f, params,Us);
R = RunGLESimulation(rp,x_0, b,c,d,e, params, Us);


%Do the continuation if subcritical is enabled and we are not doing a
%subsequent run (i.e. if we do not know the value of r at the maxxwellpt
%denoted as rc).

% Plot the solution at the end time and any vertical curves.
plotSols;
plot(x, R(end,:),'--k','linewidth',2)
plot(x, -R(end,:),'--k','linewidth',2)


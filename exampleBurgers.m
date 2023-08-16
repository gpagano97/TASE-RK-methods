%% Main code: exampleBurgers.m
Fun = @funBurgers;
Tmethod = [20 30 40]; % Select the TASE-RK methods 
nTmethods = length(Tmethod);
jacup = 0; % We want 'constant Jacobian'
if (jacup == 0)
    Jac = @jacBurgersfix;
elseif (jacup == 1)
    Jac = @jacBurgers;
end
inN = 8; finN = 12; % We do the time integration using 2^(inN),...,2^(finN) time grid points

global epsilon M L1 L2 
epsilon = 1/10; M = 32; % Spatial grid points
xspan = [0 2*pi]; tspan = [0 4];
Deltax = (xspan(2)-xspan(1))/M; % Spatial step-size

% Initial conditions
y0 = []; y0(1:M/2,1) = ones(M/2,1); y0(M/2+1:M,1) = zeros(M/2,1);

% Space-discetization: order-four FD and periodic BC
e = ones(M,1); r = 1/(12*Deltax^2);
L1 = spdiags([-r*e 16*r*e -30*r*e 16*r*e -r*e], -2:2, M, M); L1(1,M-1) = -r; L1(M-1,1) = -r; L1(M,2) = -r; L1(2,M) = -r; L1(1,end) = 16*r; L1(end,1) = 16*r;
r = 1/(12*Deltax);
L2 = spdiags([r*e -8*r*e 0*r*e 8*r*e -r*e], -2:2, M, M); L2(1,M-1) = r; L2(M-1,1) = -r; L2(M,2) = -r; L2(2,M) = r; L2(1,end) = -8*r; L2(end,1) = 8*r;

% Compute a reference solution with ode15s
options = odeset('RelTol',5e-14,'AbsTol',5e-14);
[tode15s,yode15s] = ode15s(@funBurgers,tspan,y0,options);
YrefT = yode15s(end,:)';

for nm = 1:nTmethods % We apply all the selected TASE-RK
    i = 1;
    for N = 2.^(inN:finN) % For the simulations done
        [yTTRK,yTRK,t,CPUtimeTRK] = TASERK(N,tspan,y0,Fun,Jac,Tmethod(nm),jacup); 
        errT_TRK(i,nm) = norm(yTTRK-YrefT,inf); % Error
        cd_TRK(i,nm) = -log10(errT_TRK(i,nm));
        CPUtime_TRK(i,nm) = CPUtimeTRK; % CPU time
        i = i + 1;
    end
    
    l = length(cd_TRK(:,nm));   
    for i = 2:l % Compute the estimated order
        pest_TRK(i-1,nm) = (cd_TRK(i,nm)-cd_TRK(i-1,nm))/log10(2);
    end    
end

% Print results
format short e
errT_TRK
CPUtime_TRK
format short 
pest_TRK
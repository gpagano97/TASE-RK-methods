%% Main code: exampleEuler.m
Fun = @funEuler;
Tmethod = [40]; % Select the TASE-RK method
nTmethods = length(Tmethod);
jacup = 0; % We want 'constant Jacobian'
if (jacup == 0)
    Jac = @jacEulerfix;
elseif (jacup == 1)
    Jac = @jacEuler;
end

% Initial conditions
global y0
tspan = [0 10];
y0 = [1;0;0.9];
N = 5000; % Number of grid intervals

% Compute a reference solution with ode15s
options = odeset('RelTol',5e-14,'AbsTol',5e-14);
[tode15s,yode15s] = ode15s(@funEuler,tspan,y0,options);
YrefT = yode15s(end,:)';

[yTTRK,yTRK,t,CPUtimeTRK] = TASERK(N,tspan,y0,Fun,Jac,Tmethod,jacup);

% Print results
format short e
errT_TRK = norm(yTTRK-YrefT,inf) % Error
CPUtimeTRK % CPU time
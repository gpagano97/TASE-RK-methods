function [yT,y,t,CPUtime] = TASERK(N,tspan,y0,Fun,Jac,Method,jacup)

% Fixing explicit RK tableau and TASE operator coefficients
switch Method
    case 20 % 20-midpoint method of order s=p=2
        s = 2; p = 2; alpha = [3 3/2];
        A = [0 0;1/2 0]; c = [0 1/2]; b = [0 1];

    case 30 % 30-Ralston method of order s=p=3
        s = 3; p = 3; alpha = [2.31469 1.87961 1.58222];
        A = [0 0 0;1/2 0 0;0 3/4 0];
        c = [0 1/2 3/4]; b = [2/9 1/3 4/9];

    case 40 % 40-Kutta method of order s=p=4
        s = 4; p = 4; alpha = [3.939556 2.450558 2.227083 2.061235];
        A = [0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0];
        c = [0 1/2 1/2 1]; b = [1/6 1/3 1/3 1/6];
end

% Computation of gamma
alpham1 = 1./alpha; gamma = [];
for i = 1:p
    prod = 1;
    for j = 1:p
        if i ~= j
            prod = prod*(alpham1(i)-alpham1(j));
        end
    end
    gamma(i) = alpham1(i)^(p-1)/prod;
end

% Initialization
t = linspace(tspan(1),tspan(2),N+1); % t: discrete time grid
h = (tspan(2)-tspan(1))/N; % h: constant time step-size
d = length(y0); % d: dimension of the problem
Id = eye(d); % Id: Identity matrix of order d
n = 1; % n: index representing the current step
y = y0;

if (jacup == 0) % If we want a 'constant Jacobian'
    C = cputime;
    Jn = Jac();
    for l = 1:p % Compute, outside the loop, the LU factorizations of Id-alpha(l)*(h*Jn)
        [Ll(1:d,(l-1)*d+1:l*d),Ul(1:d,(l-1)*d+1:l*d)] = lu(Id-alpha(l)*(h*Jn));
    end
    for n = 2:N+1
        Y = zeros(d,s); % Block matrix with stages in columns
        Y(:,1) = y(:,n-1);
        f = zeros(d,s); % Block matrix with h*fi in columns
        F = zeros(d,s); % Block matrix with Tp*h*fi in columns
        for i = 1:s-1 % Compute all the stages
            f(:,i) = h*Fun(t(n-1)+h*c(i),Y(:,i));
            for l = 1:p
                F(:,i) = F(:,i) + gamma(l)*(Ul(1:d,(l-1)*d+1:l*d) \ (Ll(1:d,(l-1)*d+1:l*d)\f(:,i)));
            end
            Y(:,i+1) = y(:,n-1) + sum(kron(A(i+1,:),ones(d,1)).*F,2);
        end
        f(:,s) = h*Fun(t(n-1)+h*c(s),Y(:,s));
        for l = 1:p
            F(:,s) = F(:,s) + gamma(l)*(Ul(1:d,(l-1)*d+1:l*d) \ (Ll(1:d,(l-1)*d+1:l*d)\f(:,s)));
        end
        y(:,n) = y(:,n-1) + sum(kron(b,ones(d,1)).*F,2);
    end
    Cf = cputime;
elseif (jacup == 1) % If we want exact Jacobian
    C = cputime;
    for n = 2:N+1
        Jn = Jac(t(n-1),y(:,n-1)); % Update Jn at each step
        for l = 1:p % Compute, at each step, the LU factorizations of Id-alpha(l)*(h*Jn)
            [Ll(1:d,(l-1)*d+1:l*d),Ul(1:d,(l-1)*d+1:l*d)] = lu(Id-alpha(l)*(h*Jn));
        end
        Y = zeros(d,s); % Block matrix with the stages in column
        Y(:,1) = y(:,n-1);
        f = zeros(d,s); % Block matrix with h*fi in column
        F = zeros(d,s); % Block matrix with Tp*h*fi in column
        for i = 1:s-1 % Compute all the stages
            f(:,i) = h*Fun(t(n-1)+h*c(i),Y(:,i));
            for l = 1:p
                F(:,i) = F(:,i) + gamma(l)*(Ul(1:d,(l-1)*d+1:l*d) \ (Ll(1:d,(l-1)*d+1:l*d)\f(:,i)));
            end
            Y(:,i+1) = y(:,n-1) + sum(kron(A(i+1,:),ones(d,1)).*F,2);
        end
        f(:,s) = h*Fun(t(n-1)+h*c(s),Y(:,s));
        for l = 1:p
            F(:,s) = F(:,s) + gamma(l)*(Ul(1:d,(l-1)*d+1:l*d) \ (Ll(1:d,(l-1)*d+1:l*d)\f(:,s)));
        end
        y(:,n) = y(:,n-1) + sum(kron(b,ones(d,1)).*F,2);
    end
    Cf = cputime;
end

CPUtime = Cf - C;
yT = y(:,end);
end
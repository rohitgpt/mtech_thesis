%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r1, r2 are inner and outer radius in mm
% R is the radius of curvature of the bamboo
% sigma is the max stress on any element
% M_min is the min bending moment of the bamboo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M, y_inv, min_x, max_x, x] = Prob_1(E1, E2, rho1, rho2, n, r1, r2, cur, sigma, M_min)
% E1 = 20*10^9;
% E2 = (4.25)*10^9;
% Elimit = a*E1;
% rho1 = 2000;
% rho2 = 1000;
r1=r1/1000;
% dr = @(X) (r2-r1*1000)/(n*1000);
dr = (r2-r1*1000)/(n*1000);
e = 0.001;

% R = {@(X) r1};
% R = zeros(n,1);
i = 1:n;
R = r1+dr*i;
% for i=2:n
%     R{i} = @(X) r1 + dr(X)*i;
% end

% E = @(X) R{1}(X)^3*(X(2)*E1 + X(n+2)*E2);
% rho = @(X) R{1}(X)*(X(2)*rho1+X(n+2)*rho2);
E = @(X) R(1)^3*(X(2)*E1 + X(n+2)*E2);
rho = @(X) R(1)*(X(2)*rho1+X(n+2)*rho2);
for i=2:n;
%     E = @(X) E(X) + (R{i}(X))^3*(X(i+1)*E1 + X(n+1+i)*E2);
%     rho=@(X) rho(X)+(R{i}(X))*(X(i+1)*rho1+X(n+1+i)*rho2);
    E = @(X) E(X) + (R(i))^3*(X(i+1)*E1 + X(n+1+i)*E2);
    rho=@(X) rho(X)+(R(i))*(X(i+1)*rho1+X(n+1+i)*rho2);
end

stiff = @(X) pi*E(X);
density = @(X) 2*pi*rho(X); 

% y_neg = @(X) -1*stiff(X)/density(X);
y_inv = @(X) density(X)/stiff(X);

M = {@(X) y_inv(X);};
%A=@(X) 0;
for i=2:n+1
    %r = @(X) xmin+((i-1)/n)*(1-xmin);    
%    A=@(X) A(X)+(R{i-1}(X))^3;
    M{i} = @(X) (X(i)+X(n+i)-1);
%     M{i+n} = @(X) 1*(R{i-1}(X)*(X(i)*E1 + X(n+i)*E2)/cur - sigma);
    M{i+n} = @(X) 1*(R(i-1)*(X(i)*E1 + X(n+i)*E2)/cur - sigma);
%    M{i+n} = @(X) 1*(R{i-1}(X)*(X(i)*E1 + X(n+i)*E2) - b*Elimit);
%     M{2*i+n} = @(X) 1*(R{i-1}(X)*(X(i)*E1 - X(n+i)*E2) - Elimit);
end
% M{n+2} = @(X) -1*(E(X)- cur*M_min/(dr(X)*pi));
M{n+2} = @(X) 1*(E(X)- cur*M_min/(dr*pi));
M{n+3} = @(X) y_inv(X) - 1e-3;
%M{2*n+2} = @(X) -1*(E(X) - c*A(X)*Elimit);

min_x = e*ones(1, 2*n+1);
max_x= (1-10*e)*ones(1, 2*n+1);
max_x(1) = r2/1000;
x = e*ones(1, 2*n+1);
end

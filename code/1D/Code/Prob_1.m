function [M,min_x, max_x, x] = Prob_1(n,a,b,c)
E1 = 20*10^9;
E2 = (4.25)*10^9;
Elimit = a*E1;
rho1 = 2000;
rho2 = (1.1)*1000;
e = 0.001;
xmin=0.33;
R = {@(X) xmin+(1/n)*(1-xmin);};
for i=2:n
    R{i} = @(X) xmin+(i/n)*(1-xmin);
end
%r = @(X) xmin+(1/n)*(1-xmin);   
E = @(X) R{1}(X)^3*(X(2)*E1 + X(n+2)*E2);
rho = @(X) R{1}(X)*(X(2)*rho1+X(n+2)*rho2);
for i=2:n;
    %r = @(X) xmin+(i/n)*(1-xmin);
    E = @(X) E(X) + (R{i}(X))^3*(X(i+1)*E1 + X(n+1+i)*E2);
    rho=@(X) rho(X)+(R{i}(X))*(X(i+1)*rho1+X(n+1+i)*rho2);
end

stiff = @(X) pi*E(X);
density = @(X) 2*pi*rho(X); 

% y_neg = @(X) -1*stiff(X)/density(X);
y_inv = @(X) density(X)/stiff(X);

M = {@(X) y_inv(X);};
A=@(X) 0;
for i=2:n+1
    %r = @(X) xmin+((i-1)/n)*(1-xmin);    
    A=@(X) A(X)+(R{i-1}(X))^3;
    M{i} = @(X) (X(i)+X(n+i)-1);
    M{i+n} = @(X) 1*(R{i-1}(X)*(X(i)*E1 + X(n+i)*E2) - b*Elimit);
%     M{2*i+n} = @(X) 1*(R{i-1}(X)*(X(i)*E1 - X(n+i)*E2) - Elimit);
end
M{2*n+2} = @(X) -1*(E(X) - c*A(X)*Elimit);

min_x = e*ones(1, 2*n+1);
max_x= (1-10*e)*ones(1, 2*n+1);
x = ones(1, 2*n+1);
end

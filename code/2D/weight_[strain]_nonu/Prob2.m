function [M, min_x, max_x, x] = Prob2(nx, ny, E, rho, p, max_strain)
m = length(E);
[U] = @(X) finite_element(nx, ny, E, m, X, p);
subindex = @(A, a, b) A(a, b);
[KE] = lk;
%Objective function
obj = @(X) 0;
syms xa ya;
n = 1/4*[(1-xa)*(1-ya) (1+xa)*(1-ya) (1+xa)*(1+ya) (1-xa)*(1+ya)];
B = sym('b', [3 8]);
for i=1:4
  B(1:3,i*2-1:i*2) = [diff(n(i), xa) 0; 0 diff(n(i), ya); diff(n(i), ya) diff(n(i),xa)];
end
B = subs(subs(B, xa, 0), ya, 0) + ...
subs(subs(B, xa, 1), ya, 0) + ...
subs(subs(B, xa, 0), ya, 1) + ...
subs(subs(B, xa, 1), ya, 1);
M = {@(X) 0};
for ely = 1:ny
  for elx = 1:nx
    n1 = (ny+1)*(elx-1)+ely; 
    n2 = (ny+1)* elx   +ely;
    Ue = @(X) subindex(U(X), [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
    xe = @(X) xelement(X, E, nx, ely, elx);
    e = @(X) dot(xe(X).^p, E);
    obj = @(X) obj(X) + objective(Ue(X), E, xe(X), p, rho);
    M{(ely-1)*nx+elx+1} = @(X) max(Ue(X)) - max_strain;
%     M{(ely-1)*nx+elx+1} = @(X) max(e(X)*B*Ue(X)) - max_stress;
  end
end

M{1} = @(X) obj(X);
epsilon = 0.001;
min_x = epsilon*ones(ny,m*nx);
max_x = (1-1*epsilon)*ones(ny,m*nx);
x = epsilon*ones(ny,m*nx);
end

function s = objective(Ue, E, xe, p, rho)
[KE] = lk;
s = dot(xe, rho);
% s = dot(xe.^p, rho)*dot(xe.^p, E)*(mean(xe))^(p)*Ue'*KE*Ue;
% s = dot(xe.^p, E)*Ue'*KE*Ue;
end

function xe = xelement(x, E, nx, ely, elx)
xe = [0];
for i=1:length(E)
  xe(i) = x(ely, (i-1)*nx+elx);
end
end

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=finite_element(nx,ny, E, m, x, p)
[KE] = lk; 
K = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
F = sparse(2*(ny+1)*(nx+1),1); U = zeros(2*(ny+1)*(nx+1),1);
for elx = 1:nx
  for ely = 1:ny
    n1 = (ny+1)*(elx-1)+ely; 
    n2 = (ny+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    xe = [0];
    for i=1:m
      xe(i) = x(ely, (i-1)*nx+elx);
    end
    K(edof,edof) = K(edof,edof) + dot(xe.^p, E)*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F([1+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)],1) = 0.4*1e6;
fixeddofs   = union([1:2:2*(ny+1)],union([2:(ny+1)*2:2*(nx+1)*(ny+1)],[2*(ny+1):2*(ny+1):2*(ny+1)*(nx+1)]));
alldofs     = [1:2*(ny+1)*(nx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
end
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
%Possion's ratio
nu = 0.3;
E = 1.; 
% syms nu E;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = 1/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
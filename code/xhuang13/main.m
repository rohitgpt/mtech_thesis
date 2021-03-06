function main(L, H, Vf, nx, ny)
delta = 0.1; %change in volume in every iteration
rmin = 3;     %for mesh independent filter

% x = random_init(nx, ny, Vf);   %initial design for microscale FEM
x = design1(nx, ny);
% x = design2(nx, ny);
E = [1 0.1];
syms xa ya;
alpha = 25e-6;
beta = 25e-6;
global gauss;
gauss = 0.57735;
% n = 1/(4*alpha*beta)*[(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
n = 1/(4*alpha*beta)*[(alpha-xa)*(beta+ya) (alpha+xa)*(beta+ya) ...
                      (alpha+xa)*(beta-ya) (alpha-xa)*(beta-ya)];
B = sym('b', [3 8]);
for i=1:4
  B(1:3,i*2-1:i*2) = [diff(n(i), xa)  0;...
                      0               diff(n(i), ya);...
                      diff(n(i), ya)  diff(n(i),xa);...
                      ];
end
nu=0.3;
Ddash = 1/(1-nu^2)*[1   nu  0;
                    nu  1   0;
                    0   0   (1-nu)/2;];
KE = double(int(int(B'*Ddash*B, xa, -alpha, alpha), ya, -beta, beta));
b = double(int(int(B, xa, -alpha, alpha), ya, -beta, beta));
% b = 1/4*(subs(subs(B, xa, -alpha*gauss), ya, -beta*gauss) + ...
% subs(subs(B, xa, alpha*gauss), ya, -beta*gauss) + ...
% subs(subs(B, xa, alpha*gauss), ya, beta*gauss) + ...
% subs(subs(B, xa, -alpha*gauss), ya, beta*gauss));

prev_C=1e10;
prev_s = zeros(ny, nx);
C = 1e5;
count=0;
v = sum(x(:))/(nx*ny);
while (prev_C-C)/C>0.05 || v~=Vf
  prev_C = C;
  count=count+1;
  plot_fig(x);
  u = microFEM(nx, ny, x, KE, b, Ddash, E);

  Dh = homogenization(nx, ny, x, b, u, Ddash, E);
  [C, U] = macroFEM(L, H, Dh, B);
%   disp(abs((prev_C-C)/C));
  disp(C);
%   disp(prev_C);
  s = sensitivity(nx, ny, x, L, H, Ddash, E, u, B, U);
  disp(s);
  s_filtered = apply_filter(nx, ny, s, rmin);
%   disp(s_filtered);
  avg_s = 0.5*(s_filtered + prev_s);
  prev_s = avg_s;

  if v>Vf
    if v*(1-delta)<Vf
      v = Vf;
    else
      v = v*(1-delta);
    end
  else
    if v*(1+delta)>Vf
      v = Vf;
    else
      v = v*(1+delta);
    end
  end

  [~, index] = sort(avg_s(:), 'descend');
  y = zeros(nx*ny, 1);
  y(index(1:ceil(v*(nx*ny))))=1;
  x = reshape(y, ny, nx);
end

end

function s_filtered = apply_filter(nx, ny, s, rmin)
s_filtered = zeros(size(s));
for i=1:nx
  for j=1:ny
    total = 0;
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),ny)
        weight = rmin-sqrt((i-k)^2+(j-l)^2);
        s_filtered(j, i) = s_filtered(j, i)+s(l, k)*weight;
        total = total + weight;
      end
    end
    s_filtered(j, i) = s_filtered(j, i)/total;
  end
end
end

function s = sensitivity(nx, ny, x, L, H, D, E, u, B1, U)
alpha = 1e-2; beta = 1e-2;
syms xa ya;
n = 1/(4*alpha*beta)*[(alpha-xa)*(beta+ya) (alpha+xa)*(beta+ya) ...
                      (alpha+xa)*(beta-ya) (alpha-xa)*(beta-ya)];
B = sym('b', [3 8]);
for i=1:4
  B(1:3,i*2-1:i*2) = [diff(n(i), xa)  0;...
                      0               diff(n(i), ya);...
                      diff(n(i), ya)  diff(n(i),xa);...
                      ];
end
s(1:ny, 1:nx) = 1;
for i=1:nx
  for j=1:ny
    alpha = 25e-6;    beta = 25e-6;
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    t2 = (E(1)-E(2))*double(int(int((eye(3)-B1*u(dof,:))'*D*(eye(3)-B1*u(dof,:)),...
                          xa, -alpha, alpha), ya, -beta, beta));
    alpha = 1e-2;    beta = 1e-2;
    t1 = double(int(int(B'*t2*B, xa, -alpha, alpha), ya, -beta, beta));
    temp = 0;
    for l=1:L
      for h=1:H
        n1 = (H+1)*(l-1)+h;
        n2 = (H+1)*l+h;
        dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        temp = temp + U(dof, 1)'*t1*U(dof,1);
      end
    end
    s(j, i) = x(j, i)*temp;
  end
end
end
%% Plotting
function plot_fig(x)
% disp(x);
figure(1)
colormap(gray);
imagesc(1-x);
pause(1);
end
%% Perform macro FEM 
function [C, U] = macroFEM(L, H, Dh, B)
U = sparse(2*(L+1)*(H+1), 1);
K = sparse(2*(L+1)*(H+1),2*(L+1)*(H+1));
F = sparse(2*(L+1)*(H+1), 1);
alpha = 1e-2;
beta = 1e-2;
syms xa ya;
n = 1/(4*alpha*beta)*[(alpha-xa)*(beta+ya) (alpha+xa)*(beta+ya) ...
                      (alpha+xa)*(beta-ya) (alpha-xa)*(beta-ya)];
B = sym('b', [3 8]);
for i=1:4
  B(1:3,i*2-1:i*2) = [diff(n(i), xa)  0;...
                      0               diff(n(i), ya);...
                      diff(n(i), ya)  diff(n(i),xa);...
                      ];
end
% disp(B);
%KE is calculated every iteration as Dh is changes every iteration
KE = double(int(int(B'*Dh*B, xa, -alpha, alpha), ya, -beta, beta));
% disp(KE);
for i=1:L
  for j=1:H
    n1 = (H+1)*(i-1)+j;
    n2 = (H+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    %Every element is identical, KE is assembled directly;
    K(dof, dof) = K(dof, dof) + KE;
  end
end
n2 = (H+1)*(L)+ceil(H/2);
F(2*n2+2, 1) = -1;
fixed_dofs = 1:2*(H+1);
free_dofs = setdiff(1:2*(H+1)*(L+1), fixed_dofs);
U(free_dofs, :) = K(free_dofs, free_dofs)\F(free_dofs, :);
U(fixed_dofs, :) = 0;
% disp(det(K(free_dofs, free_dofs)));
% disp(U);
C = double(0.5*F'*U);
end
%% Return homogenized elasticity matrix
function Dh = homogenization(nx, ny, x, b, u, D, E)
Dh =  zeros(3,3); %initialising Dh to be 3x3
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    Dh = Dh + (x(j, i)*E(1)+(1-x(j, i))*E(2))*(eye(3)-b*u(dof, :));         %
  end
end
Dh = double(Dh*D/(nx*ny));
end
%% Perform microFEM
function u = microFEM(nx, ny, x, ke, b, D, E)
le = b'*D;
%sparse matrices are causing error
load = zeros(2*(nx+1)*(ny+1),3);
k = zeros(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));

% load = sparse(2*(nx+1)*(ny+1), 3);
% k = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
% u = sparse(2*(nx+1)*(ny+1), 3);

for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
%     disp(dof);
    %k is the assembly of elemental stiffness matrices
    %every elemental stiffness matrix is identical, i.e. ke is calculated
    %once and for different element it is multiplied by the relative
    %density of element
    k(dof, dof) = k(dof, dof) + (x(j, i)*E(1)+(1-x(j, i))*E(2))*ke;
%     det(k(dof,dof))
    %load is assembly of elemental loads corresponding to the eqn (5) in
    %the paper. le (elemental load) is calculated once and multiplied by
    %relative density for every element.
    load(dof, :) = load(dof, :) + (x(j, i)*E(1)+(1-x(j, i))*E(2))*le;
  end
end
% disp(det(ke));
% disp(det(k));
u = k\load;
% disp(u);
end
%% Initial designs for the micro FEM
function x = random_init(nx, ny, Vf)
xinit = rand(ny, nx);
[~, index] = sort(xinit(:), 'descend');
v = Vf*(nx*ny);
y = zeros(nx*ny, 1);
y(index(1:ceil(v)))=1;
x = reshape(y, ny, nx);
end

function xinit = design1(nx, ny)
xinit(1:ny, 1:nx) = 1;
xmid = ceil(nx/2); ymid=ceil(ny/2);
xinit(ymid:ymid+1,xmid:xmid+1) = 0;
end

function xinit = design2(nx, ny)
xinit(1:ny, 1:nx) = 1;
xinit(1,1)=0;
xinit(ny,1)=0;
xinit(1,nx)=0;
xinit(ny, nx)=0;
end
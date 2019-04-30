function main(L, H, Vf, nx, ny)

delta = 0.02; %change in volume in every iteration
rmin = 4;     %for mesh independent filter

x = random_init(nx, ny, Vf);   %initial design for microscale FEM
% xinit = design1(nx, ny);
% xinit = design2(nx, ny);

syms xa ya;
alpha = 25e-6;
beta = 25e-6;
global gauss;
gauss = 0.57735;
n = 1/(4*alpha*beta)*[(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
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
b = 1/4*(subs(subs(B, xa, -alpha*gauss), ya, -beta*gauss) + ...
subs(subs(B, xa, alpha*gauss), ya, -beta*gauss) + ...
subs(subs(B, xa, alpha*gauss), ya, beta*gauss) + ...
subs(subs(B, xa, -alpha*gauss), ya, beta*gauss));

prev_C=1e10;
C = 1e5;
while (prev_C-C)/C>0.05
plot_fig(x);
u = microFEM(nx, ny, x, KE, b, Ddash);

Dh = homogenization(nx, ny, b, u, D);

[C, U] = macroFEM(L, H, Dh, B);

s = sensitivity(L, H, Ddash, b, u, B, U);

s_filtered = apply_filter(s, rmin);

avg_s = 0.5*(s_filtered + prev_s);
prev_s = avg_s;

v = sum(x(:))/(nx*ny);
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
v = v*(nx*ny);
y = zeros(nx*ny, 1);
y(index(1:v))=1;
x = reshape(y, nx, ny);

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
        s_filtered(j, i) = s(l, k)*weight;
        total = total + weight;
      end
    end
    s_filtered(j, i) = s_filtered(j, i)/total;
  end
end
end

function s = sensitivity(L, H, D1, b, u, B, U)
D2 = zeros(3,3);
s(1:ny, 1:nx) = 1;
for i=1:nx
  for j=1:ny
    alpha = 25e-6;
    beta = 25e-6;
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    s(j, i) = double(int(int((eye(3)-B*u(dof,3))'*(D1-D2)*(eye(3)-B*u(dof,3)), xa, -alpha, alpha), ya, -beta, beta));
    temp = 0;
    for l=1:L
      for h=1:H
        alpha = 1;
        beta = 1;
        n1 = (H+1)*(l-1)+h;
        n2 = (H+1)*l+h;
        dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        t1 = double(int(int(B'*s(j,i)*B, xa, -alpha, alpha), ya, -beta, beta));
        temp = temp + U(dof, 1)'*t1*U(dof,1);
      end
    end
    s(j, i) = temp;
  end
end

for i=1:L
  for j=1:H
    
  end
end

end

function plot_fig(x)
figure(1)
colormap(gray);
imagesc(1-x);
end

function [C, U] = macroFEM(L, H, Dh, B)
U = sparse(2*(L+1)*(H+1), 1);
K = sparse(2*(L+1)*(H+1),2*(L+1)*(H+1));
F = sparse(2*(L+1)*(H+1), 1);
alpha = 1;
beta = 1;
KE = double(int(int(B'*Dh*B, xa, -alpha, alpha), ya, -beta, beta));

for i=1:L
  for j=1:H
    n1 = (H+1)*(i-1)+j;
    n2 = (H+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(dof, dof) = K(dof, dof) + KE;
  end
end
n2 = H*(L-1)+ceil(H/2);
F(2*n2+2, 1) = -1;
fixed_dofs = [1:2*(H+1)];
free_dofs = setdiff([1:2*(H+1)*(L+1)], fixed_dofs);
U(free_dofs, :) = K(free_dofs, free_dofs)\F(free_dofs, :);
U(fixed_dofs, :) = 0;

C = 0.5*F'*U;
end

function Dh = homogenization(nx, ny, b, u, D)
Dh =  zeros(3,3); %initialising Dh to be 3x3
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    Dh = Dh + (eye(3)-b*u(dof, 3);
  end
end

Dh = Dh*D/(nx*ny);
end

function u = microFEM(nx, ny, x, ke, b, D)
le = b'*D;
load = sparse(2*(nx+1)*(ny+1), 3);
k = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
u = sparse(2*(nx+1)*(ny+1), 3);
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    k(dof, dof) = k(dof, dof) + x(j, i)*ke;
    load(dof, 3) = load(dof, 3) + x(j, i)*le;
  end
end
u = k\load;
end

function xinit = random_init(nx, ny, Vf)
xinit = rand(ny, nx);
for j=1:ny
  for i=1:nx
    if(xinit(j,i)>Vf) xinit(j,i)=1;
    else xinit(j,i)=0;
    end
  end
end
end

function xinit = design1(nx, ny)
xinit(1:ny, 1:nx) = 0;
xmid = ceil(nx/2); ymid=ceil(ny/2);
xinit(ymid:ymid+1,xmid:xmid+1) = 1;
end

function xinit = design2(nx, ny)
xinit(1:ny, 1:nx) = 0;
xinit(1,1)=xinit(ny,1)=xinit(1,nx)=xinit(ny, nx)=1;
end
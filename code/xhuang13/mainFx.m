function mainFx(L, H, Vf, nx, ny)
delta = 0.1; %change in volume in every iteration
rmin = 3;     %for mesh independent filter

% x = random_init(nx, ny, Vf);   %initial design for microscale FEM
x = design1(nx, ny);
% x = design2(nx, ny);
E = [1 0.01];

a = 25e-6; b = 25e-6;
% n = {@(x,y)(a-x)*(b+y)/(4*a*b);
%       @(x,y)(a+x)*(b+y)/(4*a*b);
%       @(x,y)(a+x)*(b-y)/(4*a*b);
%       @(x,y)(a-x)*(b-y)/(4*a*b)};
B = @(x,y,a,b)[ - b - y,       0, b + y,     0,   b - y,       0, y - b,     0;
              0,   a - x,     0, a + x,       0, - a - x,     0, x - a;
              a - x, - b - y, a + x, b + y, - a - x,   b - y, x - a, y - b]/(4*a*b);
nu=0.3;
D = 1/(1-nu^2)*[1   nu  0;
                nu  1   0;
                0   0   (1-nu)/2;];
ke = @(x,y) B(x,y,a,b)'*D*B(x,y,a,b);
vx = @(a,b)[-a  b;
             a  b;
             a -b;
            -a -b];
b1 = integQuad(@(x,y)B(x,y,a,b), vx(a,b));
KE = integQuad(ke, vx(a,b));

prev_C=1e10;
prev_s = zeros(ny, nx);
prev_s1 = zeros(ny, nx);
C = 1e5;
count=0;
v = sum(x(:))/(nx*ny);
while (prev_C-C)/C>0.01 || v~=Vf
  tic;
  prev_C = C;
  count=count+1;
  u = microFEM(nx, ny, x, KE, E);

  Dh = homogenization(nx, ny, x, b1, u, D, E);
  [C, U] = macroFEM(L, H, Dh, B, vx);
%   disp(abs((prev_C-C)/C));
  disp(C);
%   disp(prev_C);
  s = sensitivity(nx, ny, x, L, H, D, E, u, B, vx, U);
%   disp(s);
  s_filtered = apply_filter(nx, ny, s, rmin);
%   disp(s_filtered);
  avg_s = 0.33*(s_filtered + prev_s+prev_s1);
  prev_s = avg_s;
  prev_s1 = prev_s;
  
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
  disp(v); 
%   disp(avg_s);
  [~, index] = sort(avg_s(:), 'descend');
  y = zeros(nx*ny, 1);
  y(index(1:ceil(v*(nx*ny))))=1;
  x = reshape(y, ny, nx);
%   disp(x);
  toc;
end
plot_fig(x);

end
%%
function valInteg = integQuad(F,vertices)
    w = [1 1 1 1];
    ptGaussRef =[-0.5774   -0.5774;
                 -0.5774    0.5774;
                  0.5774   -0.5774;
                  0.5774    0.5774];
    % Shape functions
    Psi1=@(x,y)(1-x).*(1-y)/4;
    Psi2=@(x,y)(1+x).*(1-y)/4;
    Psi3=@(x,y)(1+x).*(1+y)/4;
    Psi4=@(x,y)(1-x).*(1+y)/4;
    % Shape function derivatives
    dPsi11=@(x,y) -(1-y)/4;
    dPsi21=@(x,y) (1-y)/4;
    dPsi31=@(x,y) (1+y)/4;
    dPsi41=@(x,y) -(1+y)/4;
    dPsi12=@(x,y) -(1-x)/4;
    dPsi22=@(x,y) -(1+x)/4;
    dPsi32=@(x,y) (1+x)/4;
    dPsi42=@(x,y) (1-x)/4;
    % Gradient matrix
    Jacb =@(x,y) [dPsi11(x,y), dPsi21(x,y),dPsi31(x,y),dPsi41(x,y);...
                  dPsi12(x,y), dPsi22(x,y),dPsi32(x,y),dPsi42(x,y)];
    % evaluate Shape functions on Gaussian reference points
    xi = ptGaussRef(:,1);
    eta = ptGaussRef(:,2);
    evalPsi1 = Psi1(xi,eta);
    evalPsi2 = Psi2(xi,eta);
    evalPsi3 = Psi3(xi,eta);
    evalPsi4 = Psi4(xi,eta);
    % from the change of variables function
    ptGaussDomain = [evalPsi1,evalPsi2,evalPsi3,evalPsi4]*vertices;
    % evaluate Jacobian contribution for each point
    for i=1:size(xi,1)
        evalDetJacb(i) = abs(det(Jacb(xi(i),eta(i))*vertices));
    end
    %evaluate the function on the domain points
    %evalF=F(ptGaussDomain(:,1),ptGaussDomain(:,2));
    % Finally, apply Gauss formula
    suma=zeros(size(F(vertices(1,1), vertices(1,2))));
    for i=1:size(ptGaussDomain,1)
        suma=suma+w(i)*F(ptGaussDomain(i,1),ptGaussDomain(i,2))*evalDetJacb(i);
    end
    valInteg = suma;
end
%%
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
%%
function s = sensitivity(nx, ny, x, L, H, D, E, u, B, vx, U)                %
s(1:ny, 1:nx) = 1;
for i=1:nx
  for j=1:ny
    a = 25e-6; b = 25e-6;
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    t2 = (E(1)-E(2))*integQuad(@(x,y)(eye(3)-B(x,y,a,b)*u(dof,:))'*D*(eye(3)-B(x,y,a,b)*u(dof,:)), vx(a,b));
    a = 1e-2; b = 1e-2;
    t1 = integQuad(@(x,y)B(x,y,a,b)'*t2*B(x,y,a,b), vx(a,b));
    temp = 0;
    for l=1:L
      for h=1:H
        n1 = (H+1)*(l-1)+h;
        n2 = (H+1)*l+h;
        dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        temp = temp + U(dof, 1)'*t1*U(dof,1);
      end
    end
    s(j, i) = temp;
  end
end
end
%% Plotting
function plot_fig(x)
% disp(x);
a = figure(1);
colormap(gray);
imagesc(1-x);
saveas(a,'1');
pause(1);
end
%% Perform macro FEM 
function [C, U] = macroFEM(L, H, Dh, B, vx)
U = sparse(2*(L+1)*(H+1), 1);
K = sparse(2*(L+1)*(H+1),2*(L+1)*(H+1));
F = sparse(2*(L+1)*(H+1), 1);
a = 1e-2; b = 1e-2;
syms x y;
% disp(B(x,y,a,b));
%KE is calculated every iteration as Dh is changes every iteration
KE = integQuad(@(x,y)B(x,y,a,b)'*Dh*B(x,y,a,b), vx(a,b));
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
function Dh = homogenization(nx, ny, x, b1, u, D, E)
Dh =  zeros(3,3); %initialising Dh to be 3x3
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    Dh = Dh + (x(j, i)*E(1)+(1-x(j, i))*E(2))*(eye(3)-b1*u(dof, :));         %
  end
end
Dh = double(Dh*D/(nx*ny));
end
%% Perform microFEM
function u = microFEM(nx, ny, x, ke, E)
a = 25e-6; b = 25e-6;
k = zeros(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));

for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    %k is the assembly of elemental stiffness matrices
    %every elemental stiffness matrix is identical, i.e. ke is calculated
    %once and for different element it is multiplied by the relative
    %density of element
    k(dof, dof) = k(dof, dof) + (x(j, i)*E(1)+(1-x(j, i))*E(2))*ke;
  end
end

us = (4*(ny+1)-1):2*(ny+1):2*(ny+1)*nx;
vs = (4*(ny+1)):2*(ny+1):2*(ny+1)*nx;
ue = (2*(ny+1)*(nx+1)-2*ny+1):2:(2*(ny+1)*(nx+1)-2);
ve = (2*(ny+1)*(nx+1)-2*ny+2):2:(2*(ny+1)*(nx+1)-2);
un = (2*(ny+1)+1):2*(ny+1):2*(ny+1)*nx;
vn = (2*(ny+1)+2):2*(ny+1):2*(ny+1)*nx;
uw = 3:2:2*ny;
vw = 4:2:2*ny;
fixeddofs = [2*ny+1:2*ny+2, 2*(nx+1)*(ny+1)-1:2*(nx+1)*(ny+1),...
            2*(ny+1)*(nx)+1:2*(ny+1)*(nx)+2, 1:2];
freedofs = setdiff(1:2*(nx+1)*(ny+1), fixeddofs);

%------------strain=[1 0 0]--du/dx=1--dv/dx=0--du/dy=0--dv/dy=0------------
C = zeros(2*(nx+ny), 2*(nx+1)*(ny+1));
Q = zeros(2*(nx+ny), 1);
Q(1:length(uw), :)        =2*a;
C(1:length(uw), uw)       =1;     
C(1:length(uw), ue)       =-1;    l=length(uw);
Q(l+1:l+length(uw), :)    =0;
C(l+1:l+length(uw), uw+1) =1;     
C(l+1:l+length(uw), ue+1) =-1;    l=l+length(uw);
Q(l+1:l+length(un), :)    =0;
C(l+1:l+length(un), un)   =1;     
C(l+1:l+length(un), us)   =-1;    l=l+length(un);
Q(l+1:l+length(un), :)    =0;
C(l+1:l+length(un), un+1) =1;     
C(l+1:l+length(un), us+1) =-1;    l=l+length(un);
% Q(l+1:l+2, :)             =0;     %A, uA=0, vA=0;
% C(l+1:l+2, 2*ny+1:2*ny+2) = 1;    l=l+2;
% Q(l+1:l+2, :)             =0;     %D, uD=2*b*du/dy=0, vD=2*b*dv/dy=0
% C(l+1:l+2, 1:2)           = 1;    l=l+2;
% Q(l+1)                    =2*a;   %B, uB=2*a*du/dx=a, vB=2*a*dv/dx=0
% Q(l+2)                    =0;
% C(l+1:l+2, 2*(nx+1)*(ny+1)-1:2*(nx+1)*(ny+1)) = 1;    l=l+2;
% Q(l+1)                    =2*a;   %C, uC=2*a*du/dx+2*b*du/dy=a, vC=2*a*dv/dx+2*b*dv/dy=0
% Q(l+2)                    =0;
% C(l+1:l+2, 2*(ny+1)*(nx)+1:2*(ny+1)*(nx)+2) = 1;    l=l+2;
u1 = zeros(2*(nx+1)*(ny+1),1);
u1(fixeddofs, :) = [0,0,2*a,0,2*a,0,0,0]';
alpha=1e8;
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
u1(freedofs, :)=kdash\(C(:, freedofs)'*Q);

%------------strain=[0 1 0]--du/dx=0--dv/dx=0--du/dy=0--dv/dy=1------------
C = zeros(2*(nx+ny)-8, 2*(nx+1)*(ny+1));
Q = zeros(2*(nx+ny)-8, 1);
Q(1:length(uw), :)        =0;
C(1:length(uw), uw)       =1;     
C(1:length(uw), ue)       =-1;    l=length(uw);
Q(l+1:l+length(uw), :)    =0;
C(l+1:l+length(uw), uw+1) =1;     
C(l+1:l+length(uw), ue+1) =-1;    l=l+length(uw);
Q(l+1:l+length(un), :)    =0;
C(l+1:l+length(un), un)   =1;     
C(l+1:l+length(un), us)   =-1;    l=l+length(un);
Q(l+1:l+length(un), :)    =2*b;
C(l+1:l+length(un), un+1) =1;     
C(l+1:l+length(un), us+1) =-1;    l=l+length(un);
% Q(l+1:l+2, :)             =0;     %A, uA=0, vA=0;
% C(l+1:l+2, 2*ny+1:2*ny+2) = 1;    l=l+2;
% Q(l+1)                    =0;     %D, uD=2*b*du/dy=0, vD=2*b*dv/dy=2*b
% Q(l+2)                    =2*b;   
% C(l+1:l+2, 1:2)           = 1;    l=l+2;
% Q(l+1)                    =0;     %B, uB=2*a*du/dx=0, vB=2*a*dv/dx=0
% Q(l+2)                    =0;
% C(l+1:l+2, 2*(nx+1)*(ny+1)-1:2*(nx+1)*(ny+1)) = 1;    l=l+2;
% Q(l+1)                    =0;   %C, uC=2*a*du/dx+2*b*du/dy=0, vC=2*a*dv/dx+2*b*dv/dy=2*b
% Q(l+2)                    =2*b;
% C(l+1:l+2, 2*(ny+1)*(nx)+1:2*(ny+1)*(nx)+2) = 1;    l=l+2;
u2 = zeros(2*(nx+1)*(ny+1),1);
u2(fixeddofs, :) = [0,0,0,0,0,2*b,0,2*b]';
alpha=1e8;
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
u2(freedofs, :)=kdash\(C(:, freedofs)'*Q);
%------------strain=[0 0 1]--du/dx=0--dv/dx=0.5--du/dy=0.5--dv/dy=0------------
C = zeros(2*(nx+ny), 2*(nx+1)*(ny+1));
Q = zeros(2*(nx+ny), 1);
Q(1:length(uw), :)        =0;
C(1:length(uw), uw)       =1;     
C(1:length(uw), ue)       =-1;    l=length(uw);
Q(l+1:l+length(uw), :)    =a;
C(l+1:l+length(uw), uw+1) =1;     
C(l+1:l+length(uw), ue+1) =-1;    l=l+length(uw);
Q(l+1:l+length(un), :)    =b;
C(l+1:l+length(un), un)   =1;     
C(l+1:l+length(un), us)   =-1;    l=l+length(un);
Q(l+1:l+length(un), :)    =2*b;
C(l+1:l+length(un), un+1) =1;     
C(l+1:l+length(un), us+1) =-1;    l=l+length(un);
% Q(l+1:l+2, :)             =0;     %A, uA=0, vA=0;
% C(l+1:l+2, 2*ny+1:2*ny+2) = 1;    l=l+2;
% Q(l+1)                    =b;     %D, uD=2*b*du/dy=b, vD=2*b*dv/dy=0
% Q(l+2)                    =0;   
% C(l+1:l+2, 1:2)           = 1;    l=l+2;
% Q(l+1)                    =0;     %B, uB=2*a*du/dx=0, vB=2*a*dv/dx=a
% Q(l+2)                    =a;
% C(l+1:l+2, 2*(nx+1)*(ny+1)-1:2*(nx+1)*(ny+1)) = 1;    l=l+2;
% Q(l+1)                    =b;     %C, uC=2*a*du/dx+2*b*du/dy=b, vC=2*a*dv/dx+2*b*dv/dy=a
% Q(l+2)                    =a;
% C(l+1:l+2, 2*(ny+1)*(nx)+1:2*(ny+1)*(nx)+2) = 1;    l=l+2;
u3 = zeros(2*(nx+1)*(ny+1),1);
u3(fixeddofs, :) = [0,0,0,a,b,a,b,0]';
alpha=1e8;
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
u3(freedofs, :)=kdash\(C(:, freedofs)'*Q);

u = [u1 u2 u3];
end
% %% Perform microFEM
% function u = microFEM(nx, ny, x, ke, b1, D, E)
% le = b1'*D;
% %sparse matrices are causing error
% load = zeros(2*(nx+1)*(ny+1),3);
% k = zeros(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
% 
% % load = sparse(2*(nx+1)*(ny+1), 3);
% % k = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
% % u = sparse(2*(nx+1)*(ny+1), 3);
% 
% for i=1:nx
%   for j=1:ny
%     n1 = (ny+1)*(i-1)+j;
%     n2 = (ny+1)*i+j;
%     dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
% %     disp(dof);
%     %k is the assembly of elemental stiffness matrices
%     %every elemental stiffness matrix is identical, i.e. ke is calculated
%     %once and for different element it is multiplied by the relative
%     %density of element
%     k(dof, dof) = k(dof, dof) + (x(j, i)*E(1)+(1-x(j, i))*E(2))*ke;
% %     det(k(dof,dof))
%     %load is assembly of elemental loads corresponding to the eqn (5) in
%     %the paper. le (elemental load) is calculated once and multiplied by
%     %relative density for every element.
%     load(dof, :) = load(dof, :) + (x(j, i)*E(1)+(1-x(j, i))*E(2))*le;
%   end
% end
% % disp(det(ke));
% % disp(det(k));
% u = double(k)\double(load);
% % disp(u);
% end
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
xinit(1:ny, 1:nx) = 0;
xmid = ceil(nx/2); ymid=ceil(ny/2);
xinit(ymid:ymid+1,xmid:xmid+1) = 1;
end

function xinit = design2(nx, ny)
xinit(1:ny, 1:nx) = 1;
xinit(1,1)=0;
xinit(ny,1)=0;
xinit(1,nx)=0;
xinit(ny, nx)=0;
end
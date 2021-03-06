function mainFx(L, H, Vf, nx, ny)
tic;
E = [1 1e-4];
delta = 0.02; %change in volume in every iteration
rmin = min([7,ceil(nx/4),ceil(ny/4)]);     %for mesh independent filter
% a = 1e-2; b = 1e-2; %a, b for the shape function of elements in the macro FEM
a = 1; b = 1; %a, b for the shape function of elements in the macro FEM
am = a/nx; bm = b/ny; %size for the elements in the micro FEM is determined from using a,b,nx,ny
% x = random_init(nx, ny, Vf);   %initial design for microscale FEM
% x = design1(nx, ny);
x = design2(nx, ny);
% x = design3(nx, ny);

B = @(x,y,a,b)[ y - b,     0,   b - y,       0, b + y,     0, - b - y,       0;
               0, x - a,       0, - a - x,     0, a + x,       0,   a - x;
               x - a, y - b, - a - x,   b - y, a + x, b + y,   a - x, - b - y]/(4*a*b);

nu=0.3;
D = 1/(1-nu^2)*[1   nu  0;
                nu  1   0;
                0   0   (1-nu)/2;];
KE = @(x,y) B(x,y,am,bm)'*D*B(x,y,am,bm);
vx = @(a,b)[-a  -b;
             a  -b;
             a   b;
            -a   b];
b1 = integQuad(@(x,y)B(x,y,am,bm), vx(am,bm));
ke = integQuad(KE, vx(am,bm));

prev_C=1e2;
nsteps = 2;
prev_s = zeros(ny, nx, nsteps);
C = 1e5;
count=0;
f = 1;
v = sum(x(:))/(nx*ny);
Cvec = 0;
vvec = 0;
while f>0.1 || v~=Vf || count<nsteps
  prev_C = C;
  count=count+1;
  u = microFEM(nx, ny, am, bm, x, ke, b1, D, E);

  Dh = homogenization(nx, ny, x, am, bm, b1, u, D, E, B, vx);
  [C, U] = macroFEM(L, H, a, b, Dh, B, vx);
  Cvec = [Cvec, C];
  vvec = [vvec, v];  
  s = sensitivity(nx, ny, a, b, x, L, H, D, E, u, B, vx, U);
  s_filtered = apply_filter(nx, ny, s, rmin);
  avg_s = (s_filtered + sum(prev_s,3))/(1+nsteps);
  prev_s(:,:,1) = avg_s;
  for i=2:nsteps
    prev_s(:,:,i) = prev_s(:,:,i-1);
  end
  f = abs(C-prev_C)/C*100;
%   f = abs(sum(avg_s(:)) - sum(sum(prev_s(:,:,nsteps))))
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
  toc;
  plot_fig(x, count);
end
save('Dh.mat', 'Dh');
plot_asmb(x);
% Cvec
figure
hold on;
yyaxis left
plot(1:length(Cvec)-1, full(Cvec(2:length(Cvec))), '-ok', 'LineWidth', 2, 'MarkerSize', 6); 
xlabel('Iteration');
ylabel('C, mean compliance');
yyaxis right
plot(1:length(vvec)-1, vvec(2:length(vvec)), '-sk', 'LineWidth', 2, 'MarkerSize', 6);
ylabel('Volume fraction');
legend('Compliance', 'Volume Fraction');
savefig('history');
hold off;
end

function plot_asmb(x)
n = 3;
s = size(x);
t = zeros(3*s);
for i=1:n
  for j=1:n
    t(1+(i-1)*s(1):i*s(1), 1+(j-1)*s(2):j*s(2)) = x;
  end
end
plot_fig(t,0);
end

function plot_result(nx, ny, x, u)
[X, Y] = meshgrid(1:nx+1, 1:ny+1);
Y = flipud(Y);
X1 = X + 1*reshape(u(1:2:size(u)), ny+1, nx+1);
Y1 = Y + 1*reshape(u(2:2:size(u)), ny+1, nx+1);
% disp(u);
C = ones(ny+1, nx+1);
C(1:ny, 1:nx) = 1-x;
figure(2);
pcolor(X1, Y1, C);
colormap(gray);
pause(.1);
% plot(X,Y);
end

function plot_macro(nx, ny, x, u)
[X, Y] = meshgrid(1:nx+1, 1:ny+1);
Y = flipud(Y);
X1 = X + 1*reshape(u(1:2:size(u)), ny+1, nx+1);
Y1 = Y + 1*reshape(u(2:2:size(u)), ny+1, nx+1);
% disp(u);
C = ones(ny+1, nx+1);
C(1:ny, 1:nx) = 1-x;
figure(2);
pcolor(X1, Y1, C);
colormap(gray);
% plot(X,Y);
end

%%
function valInteg = integQuad(F,vertices)
    w = [1 1 1 1];
    ptGaussRef =[-0.5774   -0.5774;
                  0.5774   -0.5774;
                  0.5774    0.5774;
                 -0.5774    0.5774];
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
    for p=i-rmin:i+rmin
      for q=j-rmin:j+rmin
        if(p<1) 
          k=nx+p; 
        elseif p>nx 
          k=p-nx;
        else
          k=p;
        end
        if(q<1) 
          l=ny+q;
        elseif q>ny 
          l=q-ny;
        else
          l=q;
        end          
        weight = rmin-sqrt((i-p)^2+(j-q)^2);
        s_filtered(j, i) = s_filtered(j, i)+s(l, k)*weight;
        total = total + weight;
      end
    end
    s_filtered(j, i) = s_filtered(j, i)/total;
  end
end
end
%%
function s = sensitivity(nx, ny, a, b, x, L, H, D, E, u, B, vx, U)                %
s(1:ny, 1:nx) = 1;
am = a/nx; bm = b/ny;
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
%     dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    dof = [2*n1+1; 2*n1+2;  2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1;];
%     dof = [2*n2+1; 2*n2+2; 2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n1+1; 2*n1+2;];
    t2 = (E(1)-E(2))*integQuad(@(x,y)(eye(3)-B(x,y,am,bm)*u(dof,:))'*D*(eye(3)-B(x,y,am,bm)*u(dof,:)), vx(am,bm));
%     t2 = D*(E(1)-E(2))*integQuad(@(x,y) (eye(3)-B(x,y,am,bm)*u(dof,:)), vx(am,bm));
    t1 = integQuad(@(x,y) B(x,y,a,b)'*t2*B(x,y,a,b), vx(a,b));
    temp = 0;
    for l=1:L
      for h=1:H
        n1 = (H+1)*(l-1)+h;
        n2 = (H+1)*l+h;
%         dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        dof = [2*n1+1; 2*n1+2;  2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1;];
%         dof = [2*n2+1; 2*n2+2; 2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n1+1; 2*n1+2;];
        temp = temp + U(dof, 1)'*t1*U(dof,1);
      end
    end
    s(j, i) = x(j,i)*temp;
  end
end
end
%% Plotting
function plot_fig(x, i)
% a = figure(1);
a = figure('visible', 'off');
colormap(gray);
imagesc(1-x);
saveas(a,[num2str(i), '.jpg']);
if i==0
    savefig(num2str(i));
end
end
%% Perform macro FEM 
function [C, U] = macroFEM(L, H, a, b, Dh, B, vx)
U = sparse(2*(L+1)*(H+1), 1);
K = sparse(2*(L+1)*(H+1),2*(L+1)*(H+1));
F = sparse(2*(L+1)*(H+1), 1);

%KE is calculated every iteration as Dh is changes every iteration
KE = integQuad(@(x,y)B(x,y,a,b)'*Dh*B(x,y,a,b), vx(a,b));

for i=1:L
  for j=1:H
    n1 = (H+1)*(i-1)+j;
    n2 = (H+1)*i+j;
%     dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    dof = [2*n1+1; 2*n1+2;  2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1;];
    %Every element is identical, KE is assembled directly;
    K(dof, dof) = K(dof, dof) + KE;
  end
end
[F, fixed_dofs] = cantilever(L, H);
% [F, fixed_dofs] = mbb(L, H);
free_dofs = setdiff(1:2*(H+1)*(L+1), fixed_dofs);
U(free_dofs, :) = K(free_dofs, free_dofs)\F(free_dofs, :);
U(fixed_dofs, :) = 0;
% plot_macro(L, H, 0, U);
C = double(0.5*F'*U);
end
%% Return homogenized elasticity matrix
function Dh = homogenization(nx, ny, x, am, bm, b1, u, D, E, B, vx)
Dh =  zeros(3,3); %initialising Dh to be 3x3
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
%     dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    dof = [2*n1+1; 2*n1+2;  2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1;];
%     dof = [2*n2+1; 2*n2+2; 2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n1+1; 2*n1+2;];
%     Dh = Dh + (x(j, i)*E(1)+(1-x(j, i))*E(2))*(eye(3)-b1*u(dof, :));         %
    Dh = Dh + (x(j, i)*E(1)+(1-x(j, i))*E(2))*integQuad(@(x,y)(eye(3)-B(x,y,am,bm)*u(dof,:))'*D*(eye(3)-B(x,y,am,bm)*u(dof,:)), vx(am,bm));
%     Dh = Dh + (x(j, i)*E(1)+(1-x(j, i))*E(2))*integQuad(@(x,y) D*(eye(3)-B(x,y,am,bm)*u(dof,:)), vx(am,bm));
  end
end
Dh = double(Dh/(4*am*bm*nx*ny));
% Dh = 0.5*(Dh+Dh');
Dh(3, [1 2]) = 0;
Dh([1 2],3) = 0
end
%% Perform microFEM
function u = microFEM(nx, ny, am, bm, x, ke, b1, D, E)
a = am; b = bm;
k = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
F = sparse(2*(nx+1)*(ny+1), 3);
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1+1; 2*n1+2; 2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1;];
%     dof = [2*n2+1; 2*n2+2; 2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n1+1; 2*n1+2;]
    %k is the assembly of elemental stiffness matrices
    %every elemental stiffness matrix is identical, i.e. ke is calculated
    %once and for different element it is multiplied by the relative
    %density of element
    k(dof, dof) = k(dof, dof) + (x(j, i)*E(1)+(1-x(j, i))*E(2))*ke;
    F(dof, :) = F(dof, :) + (x(j, i)*E(1)+(1-x(j, i))*E(2))*b1'*D;
%     F(dof, :) = F(dof, :) + x(j,i)*b1'*D;
  end
end

% u1=pbc(nx, ny, a,b,k,[1,0,0]);
% u2=pbc(nx, ny, a,b,k,[0,1,0]);
% u3=pbc(nx, ny, a,b,k,[0,0,1]);
u1=new_pbc(nx, ny, a,b,k,F,[1,0,0]);
u2=new_pbc(nx, ny, a,b,k,F,[0,1,0]);
u3=new_pbc(nx, ny, a,b,k,F,[0,0,1]);
u = [u1 u2 u3];
% plot_result(nx, ny, x, u(:,3));
end
%%
function u = new_pbc(nx, ny, a, b, k, F, strain)
us = (4*(ny+1)-1):2*(ny+1):2*(ny+1)*nx;
ue = (2*(ny+1)*(nx+1)-2*ny+1):2:(2*(ny+1)*(nx+1)-2);
un = (2*(ny+1)+1):2*(ny+1):2*(ny+1)*nx;
uw = 3:2:2*ny;

if strain(3)==1
fixeddofs = [2*ny+1:2*ny+2, 2*(nx+1)*(ny+1)-1:2*(nx+1)*(ny+1),...
            2*(ny+1)*(nx)+1:2*(ny+1)*(nx)+2, 1:2, uw+1, ue+1, un, us];

C = zeros((nx+ny-2), 2*(nx+1)*(ny+1));
Q = zeros((nx+ny-2), 1);
Q(1:length(uw), :)        =0;
C(sub2ind(size(C),1:length(uw), ue))       =1;    
C(sub2ind(size(C),1:length(uw), uw))       =-1;    l=length(uw);  
Q(l+1:l+length(un), :)    =0;
C(sub2ind(size(C),l+1:l+length(un), un+1))   =1;     
C(sub2ind(size(C),l+1:l+length(un), us+1))   =-1;    l=l+length(un);
else
fixeddofs = [2*ny+1:2*ny+2, 2*(nx+1)*(ny+1)-1:2*(nx+1)*(ny+1),...
            2*(ny+1)*(nx)+1:2*(ny+1)*(nx)+2, 1:2, uw, ue, un+1, us+1];

C = zeros((nx+ny-2), 2*(nx+1)*(ny+1));
Q = zeros((nx+ny-2), 1);
Q(1:length(uw), :)        =0;
C(sub2ind(size(C),1:length(uw), ue+1))       =1;    
C(sub2ind(size(C),1:length(uw), uw+1))       =-1;    l=length(uw);  
Q(l+1:l+length(un), :)    =0;
C(sub2ind(size(C),l+1:l+length(un), un))   =1;     
C(sub2ind(size(C),l+1:l+length(un), us))   =-1;    l=l+length(un);
end
freedofs = setdiff(1:2*(nx+1)*(ny+1), fixeddofs);
alpha=1e7;
r = F*strain';
u = sparse(2*(nx+1)*(ny+1),1);
u(fixeddofs, :) = 0;

% r = zeros(size(u));
% r = r - k(:, fixeddofs)*u(fixeddofs,:);
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
u(freedofs, :)=kdash\(r(freedofs,:) + C(:, freedofs)'*alpha*Q);
end
%%
function u = pbc(nx, ny, a, b, k, strain)
us = (4*(ny+1)-1):2*(ny+1):2*(ny+1)*nx;
ue = (2*(ny+1)*(nx+1)-2*ny+1):2:(2*(ny+1)*(nx+1)-2);
un = (2*(ny+1)+1):2*(ny+1):2*(ny+1)*nx;
uw = 3:2:2*ny;
fixeddofs = [2*ny+1:2*ny+2, 2*(nx+1)*(ny+1)-1:2*(nx+1)*(ny+1),...
            2*(ny+1)*(nx)+1:2*(ny+1)*(nx)+2, 1:2];
freedofs = setdiff(1:2*(nx+1)*(ny+1), fixeddofs);
alpha=1e7;

C = zeros(2*(nx+ny)-4, 2*(nx+1)*(ny+1));
Q = zeros(2*(nx+ny)-4, 1);
Q(1:length(uw), :)        =2*a*strain(1);
C(sub2ind(size(C),1:length(uw), ue))       =1;    
C(sub2ind(size(C),1:length(uw), uw))       =-1;    l=length(uw);  
Q(l+1:l+length(uw), :)    =2*a*strain(3)/2;
C(sub2ind(size(C),l+1:l+length(uw), ue+1)) =1;     
C(sub2ind(size(C),l+1:l+length(uw), uw+1)) =-1;    l=l+length(uw);
Q(l+1:l+length(un), :)    =2*b*strain(3)/2;
C(sub2ind(size(C),l+1:l+length(un), un))   =1;     
C(sub2ind(size(C),l+1:l+length(un), us))   =-1;    l=l+length(un);
Q(l+1:l+length(un), :)    =2*b*strain(2);
C(sub2ind(size(C),l+1:l+length(un), un+1)) =1;     
C(sub2ind(size(C),l+1:l+length(un), us+1)) =-1;    l=l+length(un);

u = sparse(2*(nx+1)*(ny+1),1);
u(fixeddofs,:) = [ 0, 0, ...
      2*a*strain(1), 2*a*strain(3)/2, ...
      2*a*strain(1)+2*b*strain(3)/2, 2*a*strain(3)/2+2*b*strain(2),...
      2*b*strain(3)/2, 2*b*strain(2)]'; 

r = zeros(size(u));
r = r - k(:, fixeddofs)*u(fixeddofs,:);
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
u(freedofs, :)=kdash\(r(freedofs,:) + C(:, freedofs)'*alpha*Q);
end
%% Macro FEM problem
function [F, fixed_dofs] = cantilever(L, H)
n2 = (H+1)*(L)+ceil(H/2);
F = sparse(2*(L+1)*(H+1), 1);
F(2*n2+2, 1) = -1;
fixed_dofs = 1:2*(H+1);
end

function [F, fixed_dofs] = mbb(L, H)
F = sparse(2*(L+1)*(H+1), 1);
F(2*floor((L+1)/2)*(H+1)+2, 1) = -1;
fixed_dofs = [2*H+1 2*H+2 2*(H+1)*(L+1)];
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

function xinit = design3(nx, ny)
xinit(1:ny, 1:nx) = 1;
xmid = ceil((nx+1)/2); ymid=ceil((ny+1)/2);
xinit(ymid-(ny-1)/4:ymid+(ny-1)/4,xmid-(nx-1)/4:xmid+(nx-1)/4) = 0;
end
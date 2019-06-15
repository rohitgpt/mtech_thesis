function testing2()
E = [1 0];
nx=13;
ny=13;
delta = 0.02; %change in volume in every iteration
rmin = min([7,ceil(nx/4),ceil(ny/4)]);     %for mesh independent filter
a = 1; b = 1; %a, b for the shape function of elements in the macro FEM
am = a/nx; bm = a/ny; %size for the elements in the micro FEM is determined from using a,b,nx,ny

% x = random_init(nx, ny, Vf);   %initial design for microscale FEM
% x = design1(nx, ny);
% x = design2(nx, ny);

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

x = design3(nx, ny);
% x(1:ny, 1:nx) = 1;
% x(2:4,2:4)=0;
plot_fig(x, 1);
% x(rand(), ran
tic;
u = microFEM(nx, ny, a/nx,b/ny,x,ke, b1, D, E);
toc;
Dh = homogenization(nx, ny, x, b1, u, D, E)
toc;
% plot_asmb(x);
end

function plot_fig(x, i)
% a = figure(1);
a = figure('visible', 'off');
colormap(gray);
imagesc(1-x);
saveas(a,[num2str(i), '.jpg']);
end

function plot_asmb(x)
n = 3;
s = size(x);
t = zeros(3*s);
for i=1:n
  for j=1:n
    t(1+(j-1)*s(2):j*s(2),1+(i-1)*s(1):i*s(1)) = x;
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
plot_fig(C, 5);
figure(2);
pcolor(X1, Y1, C);
colormap(gray);
pause(.1);
% plot(X,Y);
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
xinit(ymid-2-(ny-1)/4:ymid+2+(ny-1)/4,xmid-(nx-1)/4:xmid+(nx-1)/4) = 0;
end

function Dh = homogenization(nx, ny, x, b1, u, D, E)
Dh =  zeros(3,3); %initialising Dh to be 3x3
% nx = 3; ny = 3;
for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
%     dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    dof = [2*n1+1; 2*n1+2;  2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1;];
%     dof = [2*n2+1; 2*n2+2; 2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n1+1; 2*n1+2;];
%     b1*u(dof,:)
    Dh = Dh + (x(j, i)*E(1)+(1-x(j, i))*E(2))*(b1*u(dof, :));         %
  end
end
Dh = double(D*(eye(3)-Dh/4));
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
%     F(dof, :) = F(dof, :) + b1'*D;
  end
end
plot_fig(k, 2);
% plot_fig(F, 3);
z = hanging_nodes(x)
% u1=pbc(nx, ny, a,b,k,[1,0,0]);
% u2=pbc(nx, ny, a,b,k,[0,1,0]);
% u3=pbc(nx, ny, a,b,k,[0,0,1]);
u1=new_pbc(nx, ny, a,b,z,k,F,[1,0,0]);
u2=new_pbc(nx, ny, a,b,z,k,F,[0,1,0]);
u3=new_pbc(nx, ny, a,b,z,k,F,[0,0,1]);
u = [u1 u2 u3];
% plot_fig(u, 3);
plot_result(nx, ny, x, u(:,2));

end

function dof = hanging_nodes(x)
dof=[];
[ny, nx] = size(x);
x = x(1:ny,1:nx-1)+x(1:ny, 2:nx);
x = x(1:ny-1,1:nx-1)+x(2:ny,1:nx-1);
for i=1:nx-1
  for j=1:ny-1
    if ~x(j, i)
      n1 = (ny+1)*(i-1)+j;
      n2 = (ny+1)*i+j;
      dof = [dof, 2*n2+1, 2*n2+2];   
    end
  end
end
end

%%
function u = new_pbc(nx, ny, a, b, zerodofs, k, F, strain)
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
fixeddofs = [fixeddofs, zerodofs];
freedofs = setdiff(1:2*(nx+1)*(ny+1), fixeddofs);
alpha=1e2;
r = F*strain';
u = sparse(2*(nx+1)*(ny+1),1);
u(fixeddofs, :) = 0;

% r = zeros(size(u));
% r = r - k(:, fixeddofs)*u(fixeddofs,:);
plot_fig(alpha*(C(:, freedofs)'*C(:, freedofs)), 4);
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
%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
%     edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    edof = sort([2*n1+1; 2*n1+2;  2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1;]);
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
full(K)
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F([1+2*(nely+1)*nelx:2:2*(nely+1)*(nelx+1)],1) = 1;
fixeddofs   = union([1:2:2*(nely+1)],union([2:(nely+1)*2:2*(nelx+1)*(nely+1)],[2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1)]));
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
end

function [KE]=lk
E = 1.; 
nu = 0.3;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end

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
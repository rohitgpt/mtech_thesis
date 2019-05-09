function u = testing(nx, ny, am, bm, x, E)
a = am; b = bm;
B = @(x,y,a,b)[ - b - y,       0, b + y,     0,   b - y,       0, y - b,     0;
              0,   a - x,     0, a + x,       0, - a - x,     0, x - a;
              a - x, - b - y, a + x, b + y, - a - x,   b - y, x - a, y - b]/(4*a*b);
vx = @(a,b)[-a  b;
             a  b;
             a -b;
            -a -b];
nu=0.3;
D = 1/(1-nu^2)*[1   nu  0;
                nu  1   0;
                0   0   (1-nu)/2;];
KE = @(x,y) B(x,y,a,b)'*D*B(x,y,a,b);

k = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
ke = integQuad(KE, vx(a,b));

for i=1:nx
  for j=1:ny
    n1 = (ny+1)*(i-1)+j;
    n2 = (ny+1)*i+j;
    dof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2]
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
alpha=1e5;

%------------strain=[1 0 0]--du/dx=1--dv/dx=0--du/dy=0--dv/dy=0------------
C = sparse(2*(nx+ny), 2*(nx+1)*(ny+1));
Q = sparse(2*(nx+ny), 1);
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
u1 = sparse(2*(nx+1)*(ny+1),1);
u1(fixeddofs, :) = [0,0,2*a,0,2*a,0,0,0]';
r = zeros(size(u1));
r = r - k(:, fixeddofs)*u1(fixeddofs,:);
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
% disp(det(kdash));
u1(freedofs, :)=kdash\(r(freedofs,:) + C(:, freedofs)'*Q);
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
u2 = sparse(2*(nx+1)*(ny+1),1);
u2(fixeddofs, :) = [0,0,0,0,0,2*b,0,2*b]';
r = zeros(size(u2));
r = r - k(:, fixeddofs)*u2(fixeddofs,:);
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
u2(freedofs, :)=kdash\(r(freedofs,:) + C(:, freedofs)'*Q);
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
u3 = sparse(2*(nx+1)*(ny+1),1);
u3(fixeddofs, :) = [0,0,0,a,b,a,b,0]';
r = zeros(size(u3));
r = r - k(:, fixeddofs)*u3(fixeddofs,:);
kdash = k(freedofs, freedofs)+alpha*(C(:, freedofs)'*C(:, freedofs));
u3(freedofs, :)=kdash\(r(freedofs,:) + C(:, freedofs)'*Q);

u = [u1 u2 u3];

plot_result(nx, ny, x, u(:,1));
end

function plot_result(nx, ny, x, u)
[X, Y] = meshgrid(1:nx+1, 1:ny+1);
C = zeros(ny+1, nx+1);
C(1:ny, 1:nx) = x;
colormap(gray);
pcolor(X, Y, C);
plot(X,Y);
end
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

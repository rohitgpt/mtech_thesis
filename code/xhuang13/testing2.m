function testing2()
E = [1 1e-4]*10^2;
nx=1;
ny=1;
delta = 0.02; %change in volume in every iteration
rmin = min([7,ceil(nx/4),ceil(ny/4)]);     %for mesh independent filter
a = 1e-2; b = 1e-2; %a, b for the shape function of elements in the macro FEM
am = .5; bm = .5; %size for the elements in the micro FEM is determined from using a,b,nx,ny

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
b = integQuad(@(x,y)B(x,y,am,bm), vx(am,bm))
ke = integQuad(KE, vx(am,bm))
u = FE(1, 1, 1, 1);
end

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
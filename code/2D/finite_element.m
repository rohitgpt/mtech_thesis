%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
F([1+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)],1) = 1e6;
fixeddofs   = union([1:2:2*(ny+1)],union([2:(ny+1)*2:2*(nx+1)*(ny+1)],...
                    [2*(ny+1):2*(ny+1):2*(ny+1)*(nx+1)]));
alldofs     = [1:2*(ny+1)*(nx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
end
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
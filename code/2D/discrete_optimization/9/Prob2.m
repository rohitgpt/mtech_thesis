function [M, min_x, max_x] = Prob2(nx, ny, E, Edash, B, rho, max_strain, max_stress, force, KE)
subindex = @(A, a, b) A(a, b);
sindex = @(A,a) A(a);
Emat = @(X) funct(X, E, nx, ny);
rhomat = @(X) funct(X, rho, nx, ny);
[U] = @(X) finite_element(nx, ny, Emat(X), force, KE);
M = {@(X) 0};
for ely = 1:ny
  for elx = 1:nx
    n1 = (ny+1)*(elx-1)+ely; 
    n2 = (ny+1)* elx   +ely;
    Ue = @(X) subindex(U(X), [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
%     M{(ely-1)*nx+elx+1} = @(X) max(Ue(X)) - max_strain;
%     M{(ely-1)*2*nx+2*elx} = @(X) max(Ue(X)) - max_strain;
    max_principle_stress = @(X) stress_cal(index(Emat(X), elx, ely)*Edash*B*Ue(X));
    M{(ely-1)*nx+elx+1} = @(X) max_principle_stress(X) - max_stress;
    
    %     M{(ely-1)*2*nx+2*elx+1} = @(X) max(index(Emat(X), elx, ely)*Edash*B*Ue(X)) - max_stress;
  end
end
n1 = (ny+1)*(nx-1)+ny;
n2 = (ny+1)*nx+ny;
Ue = @(X) subindex(U(X), [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1)/elx*nx;
strn = @(X) B*Ue(X);
M{1} = @(X) sindex(strn(X), 1)*sum(sum(rhomat(X)));
% /sum(sum(Emat(X)));
min_x = zeros(1, nx*ny);
max_x = 2*ones(1, nx*ny);
end

function strs = stress_cal(s)
strs = (s(1)+s(2))/2+sqrt(0.25*(s(1)-s(2))^2+s(3)^2);
end

function Mat = funct(X, mat, nx, ny)
epsilon = 0.001;
Mat = epsilon*ones(ny, nx);
count = 1;
for j=1:nx
  for i=1:ny
    if X(count)>0
      Mat(i, j) = mat(X(count));
    end
    count = count+1;
  end
end
end

function s = index(mat, elx, ely)
s = mat(ely, elx);
end

function [U]=finite_element(nx,ny, Emat, force, KE)
F = sparse(2*(ny+1)*(nx+1),1); U = zeros(2*(ny+1)*(nx+1),1);
K = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
for elx = 1:nx
  for ely = 1:ny
    n1 = (ny+1)*(elx-1)+ely; 
    n2 = (ny+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
%     xe = xelement(x, E, nx, ely, elx);
%     K(edof,edof) = K(edof,edof) + dot(xe.^p, E)*KE;
    K(edof,edof) = K(edof,edof) + Emat(ely, elx)*KE;
  end
end
F([1+2*(ny+1)*nx, 2*(ny+1)*(nx+1)-1], 1) = force/(2*ny);
F([3+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)-2],1) = force/(ny);
% F([1+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)],1) = force/ny;

fixeddofs   = union([1:2:2*(ny+1)],union([2:(ny+1)*2:2*(nx+1)*(ny+1)],[2*(ny+1):2*(ny+1):2*(ny+1)*(nx+1)]));
alldofs     = [1:2*(ny+1)*(nx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
end
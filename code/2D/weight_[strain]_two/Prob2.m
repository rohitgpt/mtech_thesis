function [M, min_x, max_x, x] = Prob2(nx, ny, E, Edash, B, rho, p, max_strain, max_stress, force, KE)
m = length(E);
subindex = @(A, a, b) A(a, b);
obj = @(X) 0;
[U] = @(X) finite_element(nx, ny, E, X, p, force, KE);

M = {@(X) 0};
for ely = 1:ny
  for elx = 1:nx
    n1 = (ny+1)*(elx-1)+ely; 
    n2 = (ny+1)* elx   +ely;
    Ue = @(X) subindex(U(X), [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1)/elx*nx;
    xe = @(X) xelement(X, E, nx, ely, elx);
    e = @(X) dot(xe(X).^p, E)*Edash;
    obj = @(X) obj(X) + objective(Ue(X), E, xe(X), p, rho);
%     M{(ely-1)*nx+elx+1} =  @(X) max(abs(X(ely, elx)*E*Edash*B*Ue(X))) - max_stress;
% @(X) max(Ue(X)) - max_strain;
    M{(ely-1)*2*nx+2*elx} = @(X) max(Ue(X)) - max_strain;
    M{(ely-1)*2*nx+2*elx+1} = @(X) max(abs(X(ely, elx)*dot(xe(X).^p, E)*Edash*B*Ue(X))) - max_stress;
  end
end

M{1} = @(X) obj(X);
epsilon = 0.001;
min_x = 0.4*ones(ny,m*nx);
max_x = (1-1*epsilon)*ones(ny,m*nx);
x = 0.5*ones(ny,m*nx);
end

function s = objective(Ue, E, xe, p, rho)
s = dot(xe, rho);
end

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=finite_element(nx,ny, E, x, p, force, KE)
F = sparse(2*(ny+1)*(nx+1),1); U = zeros(2*(ny+1)*(nx+1),1);
K = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
for elx = 1:nx
  for ely = 1:ny
    n1 = (ny+1)*(elx-1)+ely; 
    n2 = (ny+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    xe = xelement(x, E, nx, ely, elx);
    K(edof,edof) = K(edof,edof) + dot(xe.^p, E)*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
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

function xe = xelement(x, E, nx, ely, elx)
xe = [0];
for i=1:length(E)
  xe(i) = x(ely, (i-1)*nx+elx);
end
end
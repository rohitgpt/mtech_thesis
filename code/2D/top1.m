%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function top1(nelx,nely,volfrac,penal,rmin);
% INITIALIZE
x(1:nely,1:2*nelx) = volfrac/2; 
loop = 0; 
change = 1.;
rho1=1200;
rho2=100;
E1 = 0.; 
E2 = 1;
change=0.2;
% START ITERATION
while change > 0.01  
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal);         
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  for ely = 1:nely
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
      x1e = x(ely, elx);
      x2e = x(ely, nelx+elx);
%       Ee =  (x1e+x2e)^penal*(x1e*E1+x2e*E2)*Ue'*KE*Ue;
      Ee =  (x1e*E1*x1e^penal+E2*x2e^penal)*Ue'*KE*Ue;
      
%Original
%       c = c + (x1e+x2e)^penal*Ue'*KE*Ue;
%       dc1(ely,elx) = -penal*(x1e+x2e)^(penal-1)*Ue'*KE*Ue;
%       dc2(ely,elx) = -penal*(x1e+x2e)^(penal-1)*Ue'*KE*Ue;
%Weight Minimization
      c = c + 1/(rho1*x1e^penal+rho2*x2e^penal);
      dc1(ely,elx) = -rho1*x1e^(penal-1)/(rho1*x1e^penal+rho2*x2e^penal)^2;
      dc2(ely,elx) = -rho2*x2e^(penal-1)/(rho1*x1e^penal+rho2*x2e^penal)^2;
%Initial THought
%       c = c + (rho1*x1e+rho2*x2e)*(x1e+x2e)^(penal-1)*Ee;
%       dc1(ely,elx) = -(x1e*rho1*penal+x2e*(rho2*(penal-1)+rho1))*(x1e+x2e)^(penal-2)*Ee;
%       dc2(ely,elx) = -(x2e*rho2*penal+x1e*(rho1*(penal-1)+rho2))*(x1e+x2e)^(penal-2)*Ee;
%Exponent raise to mixture of proportion
%       c = c + ((rho1*x1e+rho2*x2e)/(rho1+rho2))^penal*Ee;
%       dc1(ely,elx) = -penal*rho1*((x1e*rho1+x2e*rho2)/(rho1+rho2))^(penal-1)*Ee/(rho1+rho2);
%       dc2(ely,elx) = -penal*rho2*((x2e*rho2+x1e*rho1)/(rho1+rho2))^(penal-1)*Ee/(rho1+rho2);
%Exponenet raised to individual density of first and second material
%       c = c + (rho1*x1e^penal+rho2*x2e^penal)*Ee;
%       dc1(ely,elx) = -penal*rho1*x1e^(penal-1)*Ee;
%       dc2(ely,elx) = -penal*rho2*x2e^(penal-1)*Ee;
      dc = [dc1 dc2];
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc, loop, change);
% PRINT RESULTS
%   xold
%   x
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES  
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end 
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc, loop, change)  
l1 = 0; l2 = 100000; move = 0.3/loop;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*(dc./lmid).^(0.5)))));
%   for ely = 1:nely
%     for elx = 1:nelx
%       x1e = xnew(ely, elx);
%       x2e = xnew(ely, nelx+elx);
%       if (x1e+x2e)>1
%           xnew(ely, elx) = x1e/(x1e+x2e);
%           xnew(ely, nelx+elx) = x2e/(x2e+x1e);
%       end
%     end
%   end
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,2*nelx);
for j = 1:nely
  t=0;
  for i = 1:2*nelx
    if i>nelx
      t=nelx;
    end
    sum=0.0; 
    for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
      for k = max(i-floor(rmin),1+t):min(i+floor(rmin),t+nelx)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
E1 = 0.; 
E2 = 1.;
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    x1e = x(ely, elx);
    x2e = x(ely, nelx+elx);
%     K(edof,edof) = K(edof,edof) + (x1e+x2e)^penal*(x1e*E1+x2e*E2)*KE;
    K(edof,edof) = K(edof,edof) + (E1*x1e^penal+E2*x2e^penal)*KE;
  end
end
F([1+2*(nely+1)*nelx:2:2*(nely+1)*(nelx+1)],1) = 1;
fixeddofs   = union([1:2:2*(nely+1)],union([2:(nely+1)*2:2*(nelx+1)*(nely+1)],[2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1)]));
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% dispload = [1+2*(nely+1)*nelx:2:2*(nely+1)*(nelx+1)];
% fixeddofs   = union([1:2:2*(nely+1)],union([2:(nely+1)*2:2*(nelx+1)*(nely+1)],[2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1)]));
% alldofs     = [1:2*(nely+1)*(nelx+1)];
% freedofs    = setdiff(alldofs,union(fixeddofs, dispload));
% % SOLVING
% U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:)      
% U(fixeddofs,:)= 0;
% U(dispload, :) = 1;
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
nu = 0.3;
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

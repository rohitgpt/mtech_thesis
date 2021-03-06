% x          = randi([0,20], 5, 5)
% [~, index] = sort(x, 2, 'descend')                % Sorting index
% result     = zeros(size(x));                       % Create output
% % sub2ind(size(x), 1:5, index(:, 1).')
% result(sub2ind(size(x), 1:5, index(:, 1).')) = 1  % Mark largest
% % result(sub2ind(size(x), 1:5, index(:, 2).')) = 1;  % Mark 2nd largest

% x = randi(10, 3, 3)
% median(median(x))
% median(x(:))

% nx = 4;
% ny = 3;
% x = randi(10, nx, ny)
% [~, index] = sort(x(:), 'descend');
% v=5;
% y = zeros(nx*ny, 1);
% y(index(1:v))=1;
% reshape(y, nx, ny)
% y = x(:);
% y(index(v+1:length(y))) = 0;
% y = reshape(y, nx, ny);
% a = y>0

% tmp = zeros(size(x));
% for i=1:length(y)
%   if mod(i, ny+1)~=0 
%     tmp(mod(i,ny+1), (i-mod(i, ny+1))/(ny+1)+1) = y(i);
%   else
%     tmp(ny+1, (i-mod(i, ny+1))/(ny+1)) = y(i);
%   end
% end
% tmp

% dist_mat1 = randi([0 7], 3, 7)                     % Create Matrix
% [RowNrs,ColNrs] = find(dist_mat1>0 & dist_mat1<6);
% [RowNrSort, Idx] = sort(RowNrs);
% for k1 = 1:max(RowNrs)
%     Out{k1} = ColNrs(Idx(RowNrSort == k1));
% end
% ot = [Out{1}, Out{2}, Out{3}]

% syms xa ya alpha beta a b x y;
% % alpha = 25e-6;
% % beta = 25e-6;
% global gauss;
% gauss = 0.57735;
% % n = 1/(4*alpha*beta)*[(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
% n = [(alpha-xa)*(beta+ya) (alpha+xa)*(beta+ya) ...
%                       (alpha+xa)*(beta-ya) (alpha-xa)*(beta-ya)];
% B = sym('b', [3 8]);
% for i=1:4
%   B(1:3,i*2-1:i*2) = [diff(n(i), xa)  0;...
%                       0               diff(n(i), ya);...
%                       diff(n(i), ya)  diff(n(i),xa);...
%                       ];
% end
% subs(B, [alpha, beta, xa, ya], [a, b, x,y])
% nu=0.3;
% Ddash = 1/(1-nu^2)*[1   nu  0;
%                     nu  1   0;
%                     0   0   (1-nu)/2;];
% KE = double(int(int(B'*Ddash*B, xa, -alpha, alpha), ya, -beta, beta));
% b = double(int(int(B, xa, -alpha, alpha), ya, -beta, beta))
% b = subs(subs(B, xa, -alpha*gauss)+subs(B, xa, alpha*gauss), ya, -beta*gauss) + ...
%     subs(subs(B, xa, -alpha*gauss)+subs(B, xa, alpha*gauss), ya, beta*gauss)
  
% function temp()
% F=@(x,y) 2*exp(x+y.^2).*y; %present function
% xmin=0; % rectangle dimensions
% xmax=1;
% ymin=0;
% ymax=1;
% % -----------------------------------------------
% % In order to create a rectangular mesh in Matlab
% % we can make the following procedure
% % -----------------------------------------------
% 
% numDiv=4; %number of subdivisions
% hx=(xmax-xmin)/numDiv;
% hy=(ymax-ymin)/numDiv;
% xi=xmin:hx:xmax; %all points in the x-axis
% eta=ymin:hy:ymax; %all points in the y-axis
% [X,Y] = meshgrid(xi,eta); % create a grid of points
% Z = F(X,Y); % function values
% surf(X,Y,Z); % optional: is just to visualize the result
% [elem,vertex] = surf2patch(X,Y,Z); % the main variables
% numElem=size(elem,1) %total number of elements
% numVert= size(vertex,1) % total number of vertices
% 
% N=2; %number of Gauss points = NxN
% integralTot=0;
% for i=1:numElem %compute the integral on each element
%     v1=[vertex(elem(i,1),1),vertex(elem(i,1),2)];
%     v2=[vertex(elem(i,2),1),vertex(elem(i,2),2)];
%     v3=[vertex(elem(i,3),1),vertex(elem(i,3),2)];
%     v4=[vertex(elem(i,4),1),vertex(elem(i,4),2)];
%     vertices=[v1;v2;v3;v4];
%     elemInteg=integQuad(F,vertices,N);
%     integralTot=integralTot+elemInteg;
% end
% actualIntegVal= (exp(1)-1)^2 %exact known value
% errorInt=abs(actualIntegVal-integralTot) %absolute error
% end

% a =[0.4945   -0.1786   -0.3022    0.0137   -0.2472    0.1786    0.0549   -0.0137;
%    -0.1786    0.4945   -0.0137    0.0549    0.1786   -0.2472    0.0137   -0.3022;
%    -0.3022   -0.0137    0.4945    0.1786    0.0549    0.0137   -0.2472   -0.1786;
%     0.0137    0.0549    0.1786    0.4945   -0.0137   -0.3022   -0.1786   -0.2472;
%    -0.2472    0.1786    0.0549   -0.0137    0.4945   -0.1786   -0.3022    0.0137;
%     0.1786   -0.2472    0.0137   -0.3022   -0.1786    0.4945   -0.0137    0.0549;
%     0.0549    0.0137   -0.2472   -0.1786   -0.3022   -0.0137    0.4945    0.1786;
%    -0.0137   -0.3022   -0.1786   -0.2472    0.0137    0.0549    0.1786    0.4945];
%    
% b =[0.4945   -0.1786   -0.3022    0.0137   -0.2473    0.1786    0.0549   -0.0137;
%    -0.1786    0.4945   -0.0137    0.0549    0.1786   -0.2473    0.0137   -0.3022;
%    -0.3022   -0.0137    0.4945    0.1786    0.0549    0.0137   -0.2473   -0.1786;
%     0.0137    0.0549    0.1786    0.4945   -0.0137   -0.3022   -0.1786   -0.2473;
%    -0.2473    0.1786    0.0549   -0.0137    0.4945   -0.1786   -0.3022    0.0137;
%     0.1786   -0.2473    0.0137   -0.3022   -0.1786    0.4945   -0.0137    0.0549;
%     0.0549    0.0137   -0.2473   -0.1786   -0.3022   -0.0137    0.4945    0.1786;
%    -0.0137   -0.3022   -0.1786   -0.2473    0.0137    0.0549    0.1786    0.4945];
% 
% E = 1.; 
% nu = 0.3;
% k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
%    -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
% KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
%                   k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
%                   k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
%                   k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
%                   k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
%                   k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
%                   k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
%                   k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
% disp(det(a))
% disp(det(b))
% disp(det(KE))

% function temp()
% a = 25e-6; b = 25e-6;
% % n = 1/(4*alpha*beta)*[(alpha-xa)*(beta+ya) (alpha+xa)*(beta+ya) ...
% %                       (alpha+xa)*(beta-ya) (alpha-xa)*(beta-ya)];
% n = {@(x,y)(a-x)*(b+y)/(4*a*b), @(x,y)(a+x)*(b+y)/(4*a*b), @(x,y)(a+x)*(b-y)/(4*a*b), @(x,y)(a-x)*(b-y)/(4*a*b)};
% B = @(x,y)[ - b - y,       0, b + y,     0,   b - y,       0, y - b,     0;
%               0,   a - x,     0, a + x,       0, - a - x,     0, x - a;
%               a - x, - b - y, a + x, b + y, - a - x,   b - y, x - a, y - b]/(4*a*b);
% nu=0.3;
% D = 1/(1-nu^2)*[1   nu  0;
%                 nu  1   0;
%                 0   0   (1-nu)/2;];
% K = @(x,y) B(x,y)'*D*B(x,y);
% vertices = [-a  b;
%              a  b;
%              a -b;
%             -a -b];
% b = integQuad(B, vertices)
% k = integQuad(K, vertices)
% end
% 
% function valInteg = integQuad(F,vertices)
%     w = [1 1 1 1];
%     ptGaussRef =[-0.5774   -0.5774;
%                  -0.5774    0.5774;
%                   0.5774   -0.5774;
%                   0.5774    0.5774];
%     % Shape functions
%     Psi1=@(x,y)(1-x).*(1-y)/4;
%     Psi2=@(x,y)(1+x).*(1-y)/4;
%     Psi3=@(x,y)(1+x).*(1+y)/4;
%     Psi4=@(x,y)(1-x).*(1+y)/4;
%     % Shape function derivatives
%     dPsi11=@(x,y) -(1-y)/4;
%     dPsi21=@(x,y) (1-y)/4;
%     dPsi31=@(x,y) (1+y)/4;
%     dPsi41=@(x,y) -(1+y)/4;
%     dPsi12=@(x,y) -(1-x)/4;
%     dPsi22=@(x,y) -(1+x)/4;
%     dPsi32=@(x,y) (1+x)/4;
%     dPsi42=@(x,y) (1-x)/4;
%     % Gradient matrix
%     Jacb =@(x,y) [dPsi11(x,y), dPsi21(x,y),dPsi31(x,y),dPsi41(x,y);...
%                   dPsi12(x,y), dPsi22(x,y),dPsi32(x,y),dPsi42(x,y)];
%     % evaluate Shape functions on Gaussian reference points
%     xi = ptGaussRef(:,1);
%     eta = ptGaussRef(:,2);
%     evalPsi1 = Psi1(xi,eta);
%     evalPsi2 = Psi2(xi,eta);
%     evalPsi3 = Psi3(xi,eta);
%     evalPsi4 = Psi4(xi,eta);
%     % from the change of variables function
%     ptGaussDomain = [evalPsi1,evalPsi2,evalPsi3,evalPsi4]*vertices;
%     % evaluate Jacobian contribution for each point
%     for i=1:size(xi,1)
%         evalDetJacb(i) = abs(det(Jacb(xi(i),eta(i))*vertices));
%     end
%     %evaluate the function on the domain points
%     %evalF=F(ptGaussDomain(:,1),ptGaussDomain(:,2));
%     % Finally, apply Gauss formula
%     suma=zeros(size(F(vertices(1,1), vertices(1,2))));
%     for i=1:size(ptGaussDomain,1)
%         suma=suma+w(i)*F(ptGaussDomain(i,1),ptGaussDomain(i,2))*evalDetJacb(i);
%     end
%     valInteg = suma;
% end

X = [1 2 3; 1 2 3; 1 2 3];
Y = X';
mymap = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 0 0];
C = [3 4 5; 1 2 5; 5 5 5];
figure(1);
pcolor(X,Y,C);
colormap(mymap);
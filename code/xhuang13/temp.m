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

% syms xa ya;
% alpha = 25e-6;
% beta = 25e-6;
% global gauss;
% gauss = 0.57735;
% n = 1/(4*alpha*beta)*[(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
% % n = 1/(4*alpha*beta)*[(alpha-xa)*(beta+ya) (alpha+xa)*(beta+ya) ...
% %                       (alpha+xa)*(beta-ya) (alpha-xa)*(beta-ya)];
% B = sym('b', [3 8]);
% for i=1:4
%   B(1:3,i*2-1:i*2) = [diff(n(i), xa)  0;...
%                       0               diff(n(i), ya);...
%                       diff(n(i), ya)  diff(n(i),xa);...
%                       ];
% end
% nu=0.3;
% Ddash = 1/(1-nu^2)*[1   nu  0;
%                     nu  1   0;
%                     0   0   (1-nu)/2;];
% KE = double(int(int(B'*Ddash*B, xa, -alpha, alpha), ya, -beta, beta));
% b = double(int(int(B, xa, -alpha, alpha), ya, -beta, beta))
% b = subs(subs(B, xa, -alpha*gauss)+subs(B, xa, alpha*gauss), ya, -beta*gauss) + ...
%     subs(subs(B, xa, -alpha*gauss)+subs(B, xa, alpha*gauss), ya, beta*gauss)
  
  
F=@(x,y) 2*exp(x+y.^2).*y; %present function
xmin=0; % rectangle dimensions
xmax=1;
ymin=0;
ymax=1;
% -----------------------------------------------
% In order to create a rectangular mesh in Matlab
% we can make the following procedure
% -----------------------------------------------

numDiv=4; %number of subdivisions
hx=(xmax-xmin)/numDiv;
hy=(ymax-ymin)/numDiv;
xi=xmin:hx:xmax; %all points in the x-axis
eta=ymin:hy:ymax; %all points in the y-axis
[X,Y] = meshgrid(xi,eta); % create a grid of points
Z = F(X,Y); % function values
surf(X,Y,Z); % optional: is just to visualize the result
[elem,vertex] = surf2patch(X,Y,Z); % the main variables
numElem=size(elem,1) %total number of elements
numVert= size(vertex,1) % total number of vertices

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

function s_filtered = apply_filter(nx, ny, s, rmin)
s_filtered = zeros(size(s));
for i=1:nx
  for j=1:ny
    total = 0;
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),ny)
        weight = rmin-sqrt((i-k)^2+(j-l)^2);
        s_filtered(j, i) = s(l, k)*weight;
        total = total + weight;
      end
    end
    s_filtered(j, i) = s_filtered(j, i)/total;
  end
end
end
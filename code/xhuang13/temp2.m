% syms xa ya alpha beta a b x y;
syms xa ya;
alpha = 0.5;
beta = 0.5;
n = 1/(4*alpha*beta)*[(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
% n = [(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
B = sym('b', [3 8]);
for i=1:4
  B(1:3,i*2-1:i*2) = [diff(n(i), xa) 0; 0 diff(n(i), ya); diff(n(i), ya) diff(n(i),xa)];
end
nu=0.3;
Edash = 1/(1-nu^2)*[1   nu  0;
                nu  1   0;
                0   0   (1-nu)/2;];
KE = double(int(int(B'*Edash*B, xa, -alpha, alpha), ya, -beta, beta));
U = sym('u',[8 3]);
Dash = int(int(B*U, xa, -alpha, alpha), ya, -beta, beta)
u = [0 0 0;
     0 0 0;
     1 0 0;
     0 0 .5;
     1 0 .5;
     0 1 .5;
     0 0 .5;
     0 1 0;];
B;
B*u;
Dh= [ 0.2458    0.0960   -0.0048;
    0.0471    0.1569   -0.0040;
   -0.0239   -0.0345    0.0783;];
Dh=0.5*(Dh+Dh')
D = double(int(int((eye(3)-B*u), xa, -alpha, alpha), ya, -beta, beta));
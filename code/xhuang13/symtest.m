function symtest()

% K = sym('k', 5);
% u = sym('u', [5,1]);
% F = sym('f', [5,1]);
% 
% fixeddofs = [2,4];
% freedofs = setdiff(1:5, fixeddofs);
% F = F - K(:, fixeddofs)*u(fixeddofs, :);
% subs(F, ['u2', 'u4'], [0, 2])
% subs(subs(F, 'u2', 0), 'u4', 2)
% u1 = K(freedofs, freedofs)\F(freedofs,:);
% D = sym('D', 3);
% b = sym('b',[3,8],'real');
% u = sym('u',[8,3],'real');
% X = sym('X', 3,'real');
% syms D6_6;
% b(1,2:2:8) = 0;
% b(2,1:2:8) = 0;
% D(:,3) = 0;
% D(3,:) = 0;
% D(3,3) = D6_6;
% % f = b*u
% Dh = D - D*X`;
n = 2;
a = sym('a',[n 1],'real');
b = sym('b',n,'real');
p = sym('p',[n 1],'real');
q = sym('q',[n 1],'real');
d = a'*b*a+p'*b*p+q'*b*q
e = simplify(det(b*(a*a'+ p*p'+ q*q')))
subs(d, 
end
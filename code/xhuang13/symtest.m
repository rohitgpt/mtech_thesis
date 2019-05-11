function symtest()

K = sym('k', 5);
u = sym('u', [5,1]);
F = sym('f', [5,1]);

fixeddofs = [2,4];
freedofs = setdiff(1:5, fixeddofs);
F = F - K(:, fixeddofs)*u(fixeddofs, :);
subs(F, ['u2', 'u4'], [0, 2])
subs(subs(F, 'u2', 0), 'u4', 2)
u1 = K(freedofs, freedofs)\F(freedofs,:);

end
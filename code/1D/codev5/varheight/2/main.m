function main(Problem)
% opengl hardwarebasic
opengl('save', 'software');
n=30;

j=1;
m=11;
dens = zeros(n,m);
for i=1:m
  tic;
%     for j=1:10

r1 = 12;              %Inner Radius
r2 = 36;              %Outer Radius
cur = 1/(1/30*(1-(i^2-1)/m^2));           %Radius of curvature
sigma = 20e6;
% sigma = (20/i)*10^6;  %Max Stress
% sigma = (20)*(1-(i^2-1)/m^2)*10^6;  %Max Stress
M_min = 500*(1-(i^2-1)/m^2);    %Moment
E1 = 20*10^9;         %Young's Modulus of First Material
E2 = 4.25*10^9;       %Young's Modulus of Second Material
rho1 = 1000*(2);      %Density of First Material
rho2 = 1000;          %Density of Second Material
t = 2.5;

[M, ~, min_x, max_x, x] = Problem(E1, E2, rho1, rho2, n,  r1, r2, cur, sigma, M_min);

options = optimoptions('fmincon','Algorithm','interior-point');

options.MaxFunEvals = 100000;
options.MaxIter=2000;
options.OptimalityTolerance = 3e-6;
options.StepTolerance=1e-50 ;

[X, ~, flag] = fmincon(@(X) M{1}(X),x,[],[],[],[],min_x,max_x,@(X) mycon(X, M),options);

if(flag ~= 0)
    x = 1:n;
    r = r1/1000+x.*(r2/1000-r1/1000)/n;
%     r.*(X(2:n+1)*E1 + X(n+2:2*n+1)*E2)/cur;
    dr = (r2/1000-r1/1000)/n;
    f=figure('visible', 'off');
    plot(r*1000, X(2:n+1), 'b');
    axis([0 r2 0 1]);
    hold on;
    plot(r*1000, X(n+2:2*n+1),'r');
    title({['\sigma_{max} = ' num2str(sigma/1000000) ' MPa Moment = ' num2str(M_min) ' N-m  \rho_1 = ' num2str(rho1)],...
        'Normal Height'});
    xlabel('Radius (in mm)');
    ylabel('Proportions of different material');
    legend({'First Material', 'Second Material'}, 'Location', 'northwest');
    saveas(f, sprintf('a%d%d.jpg', i, j), 'jpeg');
    pause(0.1);
    hold off;

    g=figure('visible', 'off');
    ratio = (1-X(2:n+1)-X(n+2:2*n+1))./(X(n+2:2*n+1));
    plot(r*1000, 2*t*ratio);

    title(['Radius of Parenchyma cells'])
    xlabel('Culm Radius (in mm)')
    ylabel('Radius of Parenchyma cells (in um)');
    saveas(g, sprintf('b%d%d.jpg', i, j), 'jpeg');
    pause(0.1);

    h=figure('visible', 'off');
    dens(:,i) = rho1*X(2:n+1)+rho2*X(n+2:2*n+1);
    plot(r*1000, rho1*X(2:n+1)+rho2*X(n+2:2*n+1), 'g');
    title('Density variation along radius');
    xlabel('Radius (in mm)');
    ylabel('Density (kg/m^3)');
    saveas(h, sprintf('c%d.jpg',i),'jpeg');
    pause(0.1);
end
    toc;
end
    l=figure(4);
    c = gray;
    c = flipud(c);
    colormap(c);
    surf(1:m,r*1000, dens);
    title('Density variation along radius');
    ylabel('Radius (in mm)');
    zlabel('Density (kg/m^3)');
    xlabel('Height (in m)');
    savefig('Density variation with height.fig');   
end
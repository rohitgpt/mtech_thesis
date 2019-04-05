function main(Problem)
% opengl hardwarebasic
opengl('save', 'software');
n=30;

% i=10;
for i=1:10
    for j=1:10
% j = 30
%    for j=2:50
%        for k=2:10
% i=1; j=1; 
k=1;
r1 = 12;
r2 = 36;
cur = 30;
sigma = 20*10^6;
M_min = 200;
% sigma = 2*i*10^6;
% M_min = 5*j;
E1 = 20*10^9;
E2 = 4.25*10^9;
rho1 = 1000*(1+i/10);
rho2 = 1000;
t = 2.5;

    [M, y_inv, min_x, max_x, x] = Problem(E1, E2, rho1, rho2, n,  r1, r2, cur, sigma, M_min);

    options = optimoptions('fmincon','Algorithm','interior-point');
    
    options.MaxFunEvals = 100000;
    options.MaxIter=2000;
    options.OptimalityTolerance = j*1e-6;
    options.StepTolerance=1e-50 ;
    
    
    [X, ~, flag] = fmincon(@(X) M{1}(X),x,[],[],[],[],min_x,max_x,@(X) mycon(X, M),options);
%     y_inv(X)
    
    if(flag ~= 0)
%         fprintf('a=%d b=%d c=%d frac = %f\t r1 = %f\t v1 = %f\t v2 = %f\n', i/10, j/10, k/10 ,X(1),X(n+1), X(2*n+1));
        x = 1:n;
%         X(1)=;
%%%        X(2:n+1)
        r = r1/1000+x.*(r2/1000-r1/1000)/n;
        r.*(X(2:n+1)*E1 + X(n+2:2*n+1)*E2)/cur;
%%%        sum((r.^3).*(X(2:n+1)*E1 + X(n+2:2*n+1)*E2))
        dr = (r2/1000-r1/1000)/n;
%%%       cur*M_min/(dr*pi)
        f = figure('visible', 'off');
        plot(r*1000, X(2:n+1), 'b');
        axis([0 r2 0 1]);
        hold on;
        plot(r*1000, X(n+2:2*n+1),'r');
        title({['MaxStress = ' num2str(sigma/1000000) ' MPa MinMoment = ' num2str(M_min) ' N-m Radius = ' num2str(cur) ' m'],...
            'Large and Thin'});
        xlabel('Radius (in mm)');
        ylabel('Proportions of different material');
        legend({'First Material', 'Second Material'}, 'Location', 'northwest');
        saveas(f, sprintf('a%d%d.jpg', i, j), 'jpeg');
        pause(0.1);
        hold off;
        g = figure('visible', 'off');
        ratio = (1-X(2:n+1)-X(n+2:2*n+1))./(X(n+2:2*n+1));
        plot(r*1000, 2*t*ratio);
%         hold on;
%         plot(r*1000, sqrt((1+ratio)*(19^2-16.5^2)));
        title(['Radius of Parenchyma cells'])
        xlabel('Culm Radius (in mm)')
        ylabel('Radius of Parenchyma cells (in um)');
        saveas(g, sprintf('b%d%d.jpg', i, j), 'jpeg');
        pause(0.1);
%         hold off;
    end
        end
    end
end
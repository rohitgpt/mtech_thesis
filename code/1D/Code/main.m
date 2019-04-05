function main(Problem)
n=20;

for i=2:10
    for j=2:10
        for k=2:10
    [M,min_x, max_x, x] = Problem(n,i/10, j/10, k/10);

    options = optimoptions('fmincon','Algorithm','interior-point');
    options.MaxFunEvals = 10000000;
    options.MaxIter=20000;
    options.TolFun=1e-08;

    [X, ~, flag] = fmincon(@(X) M{1}(X),x,[],[],[],[],min_x,max_x,@(X) mycon(X, M),options);
    if(flag ~= 0)
        fprintf('a=%d b=%d c=%d frac = %f\t r1 = %f\t v1 = %f\t v2 = %f\n', i/10, j/10, k/10 ,X(1),X(n+1), X(2*n+1));
        x = 1/n:1/n:1;
        X(1)=0.33;
        r = X(1)+x.*(1-X(1));
        plot(r, X(2:n+1), 'b');
        axis([0 1 0 1]);
        hold on;
        plot(r, X(n+2:2*n+1),'r');
        title(['Axisymmetric Distribution of two materials, a = ' num2str(i/10) ' b = ' num2str(j/10) ' c = ' num2str(k/10)]);
        xlabel('Radius (outer radius scaled to 1)');
        ylabel('Proportions of different material');
        legend({'First Material', 'Second Material'}, 'Location', 'northwest');
        saveas(gcf, sprintf('%d%d%d.jpg', i-1,j-1,k-1), 'jpeg');
        pause(0.1);
        hold off;
    end
        end
    end
end
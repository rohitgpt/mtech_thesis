function timing()
x = 1000;
y = 500;
g = @() preAllocFcn(x,y);
h = @() zeros(x,y);
diffRunTime = timeit(g)-timeit(h)
end

function mArr = preAllocFcn(x,y)
for m = 1:x
    for n = 1:y
        mArr(m,n) = 0;
    end
end
end

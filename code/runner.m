w = 2.091;
x0 = [0,0,0,0,0,0,0,0];
xall = zeros(8);
i = 0;
while (true)
fun = @(A)root2d(w,A);
x = fsolve(fun, x0);
xall(i,:) = x;
if (x == x0)
    break;
else
x0 = x;
end

end
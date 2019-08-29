%b   b1  b2  b3  b4  ba  b1b2       b2b3        b1ba
%0   0   0   8   3   0   0  -8  -3  -8  -3   0   8   8  -8  -3
% rows: x, y, z, w, k

function A = quadtomatrix(g)
    syms b_1 b_2 b_3 b_4 b_a k x y z w
    vars = [b_1, b_2, b_3, b_4];
    coeffs = [x, y, z, w, k];
    for i = 1:5
        temp = g;  
        temp = subs(temp, {coeffs(i)}, {1});
        for j = 1:5
            if j ~= i
                temp = subs(temp, {coeffs(j)}, {0});
            end
        end
        A(i, 1) = subs(temp, {b_1, b_2, b_3, b_4, b_a}, {0, 0, 0, 0, 0});
        I = eye(5);
        for j = 1:5
            A(i, j+1) = subs(temp, {b_1, b_2, b_3, b_4, b_a}, ...
                {I(1, j), I(2, j), I(3, j), I(4, j), I(5, j)}) - A(i, 1);
        end
        for j = 1:3
            for k = (j+1):4
                A(i, 13 - (4-j)*(5-j)/2 + (k-j-1)) = subs(temp, {b_1, b_2, b_3, b_4, b_a}, ...
                    {I(1, j)+I(1, k), I(2, j)+I(2, k), I(3, j)+I(3, k), I(4, j)+I(4, k), 0})...
                    - A(i, j+1) - A(i, k+1) - A(i, 1);
            end
        end
        for j = 1:4
            A(i, 12+j) = subs(temp, {b_1, b_2, b_3, b_4, b_a}, {I(1, j), I(2, j), I(3, j), I(4, j), 1})...
                - A(i, j+1) - A(i, 6) - A(i, 1);
        end
    end
end
        
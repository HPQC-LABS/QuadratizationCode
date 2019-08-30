function g = matrixtoquad(M)
    M = int64(M);
    letters = ["x" "y" "z" "w" "k"];
    g = 'b_a*(';
    for j = 1: 5
        g = [g sprintf('%s*(0', letters(j))];
        if M(j, 6) ~= 0
            g = [g sprintf('%+d', M(j, 6))];
        end
        for k = 13:16
            if M(j, k) == 0
                continue;
            end
            if M(j, k) == 1
                g = [g sprintf('+b_%d', k-12)];
            end
            if M(j, k) == -1
                g = [g sprintf('-b_%d', k-12)];
            end
            if M(j, k) ~= 1 && M(j, k)~= -1
                g = [g sprintf('%+d*b_%d', M(j, k), k-12)];
            end
        end
        g = [g ')'];
        if j ~= 5
            g = [g '+'];
        end
    end
    g = [g sprintf(')+')];
    
    terms = [1 1 1 2 2 3;
        2 3 4 3 4 4];
    for j = 1: 5
        g = [g sprintf('%s*(0', letters(j))];
        for k = 1: 12
            if k == 6
                continue;
            end
            if M(j, k) == 0
                continue;
            end
            if M(j, k) == 1 && k ~= 1
                g = [g '+'];
            end
            if M(j, k) == -1 && k ~= 1
                g = [g '-'];
            end
            if (M(j, k) ~= 1 && M(j, k)~= -1 )||k == 1
                g = [g sprintf('%+d ', M(j, k))];
                if k ~=1
                    g = [g '*'];
                end
            end
            
            if 2<= k && k <=5
                g = [g sprintf('b_%d', k-1)];
            end
            
            if k >= 7
                g = [g sprintf('b_%d*b_%d', terms(1, k-6), terms(2, k-6))];
            end
        end
        g = [g ')'];
        if j ~= 5
            g = [g '+'];
        end
    end
    syms b_1 b_2 b_3 b_4 b_a k x y z w
    g = str2sym(g);
end
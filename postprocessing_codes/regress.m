% inputsample
%b123 b234     b1234  b   b1  b2  b3  b4  ba  b1b2       b2b3        b1ba
% -8 -3 -3 -3  6:     0   0   0   8   3   0   0  -8  -3  -8  -3   0   8   8  -8  -3
% -8 -3 -3 -3  6:     0   0   8   0   3   0  -8   0  -3  -8   0  -3   8  -8   8  -3
%

f = fopen("4varCasesNotDone.txt");

A = fscanf(f, "%f %f %f %f %f: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", [21, Inf]);
A = A';

LHS = cell(1, 8);
RHS = cell(1, 8);
quadCoeff = cell(1, 8);

for i = 1: 8
    LHS{i} = zeros(0, 5); 
    RHS{i} = zeros(0, 16);
end
m = size(A, 1);
 
for i = 1:m
    caseno = 0;
    x = A(i, 1);
    y = A(i, 2);
    z = A(i, 3); 
    w = A(i, 4);
    k = A(i, 5);
    % filter out lines of satisfying a condition, then use following:
    % explained = [explained; LHS]; explained = unique(explained, 'row')
    %

    if x<=-k && -k/2<=y && z<=0 && 0<=w && A(i, 6) == 0 && A(i, 7) == -x 
        caseno = 1;
    end

    if x<=-k && -k/2<=y && z<=0 && w<=0 &&  A(i, 7) == -A(i, 6) && A(i, 8) ~=0 && A(i, 9) ~= 0
        caseno = 2;
    end

    if -k<=x && 0<=z && A(i, 8) < 0 && A(i, 7) ~= A(i, 6) && A(i, 10)<0 
        caseno = 3;
    end

    if -k<=x && 0<=z && A(i, 6) == 0 && A(i, 7) == -x && A(i, 8) == 0
        caseno = 4;
    end
    
    if -k<=x && z<=0 && 0<=w && A(i, 6) == 0 && A(i, 7) == -x && A(i, 8) == 0
        caseno = 5;
    end
    
    if -k<=x && z<=0 && 0<=w && A(i, 7) < 0 && A(i, 6) > -A(i, 7) && A(i, 10)<0
        caseno = 6;
    end
    
    if -k<=x && w<=0 && A(i, 21) <= 0 && A(i, 20) <= 0 && A(i, 20)>= A(i, 21) && A(i, 11) > 0
        caseno = 7;
    end
    
    if -k<=x && w<=0 && A(i, 7) <= 0 && A(i, 6) == -A(i, 7) && A(i, 8) > 0 && A(i, 9) > 0
        caseno = 8; 
    end
    
    if caseno ~= 0
        LHS{caseno} = [LHS{caseno}; A(i, 1:5)];
        RHS{caseno} = [RHS{caseno}; A(i, 6:21)];
    end
end

explained = zeros(0, 5);
for i = 1:8
    quadCoeff{i} = LHS{i}\RHS{i};
    quadCoeff{i};
    explained = [explained; LHS{i}];
end

explained = unique(explained, 'row');
allcase = unique(A(:, 1:5), 'row');

isequal(explained, allcase);

fileid = fopen('temp.tex', 'w');
for caseno = 3:8
M = int64(quadCoeff{caseno});
letters = ["x" "y" "z" "w" "k"];
g = 'b_a*(';
for j = 1: 5
    g = [g sprintf('%s*(', letters(j))];
    if M(j, 6) ~= 0
        g = [g sprintf('%d', M(j, 6))];
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
    g = [g sprintf('%s*(', letters(j))];
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
syms b_1 b_2 b_3 b_4 
g
gsym = str2sym(g)

fprintf(fileid, ['\n\\begin{lemma}\n If \\textcolor{red}{???}, then $' g '$ is a quadratisation.\\\\\n \\end{lemma}\n']);

fprintf(fileid, "\n\\begin{proof}\n");
for b_all = int64(0:15)
    b(1) = mod(idivide(b_all, int64(8)), 2);
    b(2) = mod(idivide(b_all, int64(4)), 2);
    b(3) = mod(idivide(b_all, int64(2)), 2);
    b(4) = mod(b_all, 2);
    fprintf(fileid, "If $\\vc b = %d%d%d%d$, then $g = ", b(1), b(2), b(3), b(4));
    fprintf(fileid, "%s", subs(gsym, {b_1, b_2, b_3, b_4}, {b(1), b(2), b(3), b(4)}));
    fprintf(fileid, '$.\n\n');
end
fprintf(fileid, "\\end{proof}\n");
end
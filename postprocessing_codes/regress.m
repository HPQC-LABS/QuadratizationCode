% inputsample
%b123 b234     b1234  b   b1  b2  b3  b4  ba  b1b2       b2b3        b1ba
% -8 -3 -3 -3  6:     0   0   0   8   3   0   0  -8  -3  -8  -3   0   8   8  -8  -3
% -8 -3 -3 -3  6:     0   0   8   0   3   0  -8   0  -3  -8   0  -3   8  -8   8  -3
%

f = fopen("4varCasesNotDone.txt");

A = fscanf(f, "%f %f %f %f %f: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", [21, Inf]);
A = A';
totalcase = 20;
LHS = cell(1, totalcase);
RHS = cell(1, totalcase);
quadCoeff = cell(1, totalcase);
explained = zeros(0, 21);
unexplained = zeros(0, 21);


for i = 1: totalcase
    LHS{i} = zeros(0, 5);
    RHS{i} = zeros(0, 16);
end
m = size(A, 1);

isexplained = zeros(1, m);
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
    
    if A(i, 6) == 0 && A(i, 7) == -x && A(i, 8) == 0 && A(i, 9) ==0 && A(i, 12)== x ...
            && A(i, 14) == z + w + k && A(i, 10) == -y  && A(i, 20) == -x - y + z && A(i, 16) == y
        caseno = 1;
        
    end
    
    if A(i, 6) == x + y + z + w + 3*k && A(i, 7) == -x -z-w-2*k ...
             && A(i, 8) == -x - y - w -2*k && A(i, 9) == -x -y -z -2*k && A(i, 21) == y + z + w + 2*k 
        caseno = 2; %bit flipping b1, b4 from 1
        
    end

    if x <= -k && 0 <= y && A(i, 6) == 0 && A(i, 7) == 0 && A(i, 8) == 0 && A(i, 9) == 0 ...
            && A(i, 12) == w && A(i, 14) == 0
        caseno = 3; %bit flipping b2, b3, b4 from 1
    end
    
    if A(i, 6) == 0 && A(i, 7) == 0 && A(i, 8) == -x && A(i, 9) ==0 && A(i, 10)== -z ...
            && A(i, 11) == y + w + k && A(i, 12) == x && A(i, 14) == y
        caseno = 4; %bit flipping b1, b2 from 1, coeff for z wrong
        
    end
    
    
    if caseno ~= 0
        LHS{caseno} = [LHS{caseno}; A(i, 1:5)];
        RHS{caseno} = [RHS{caseno}; A(i, 6:21)];
        isexplained(i) = 1;
        explained = [explained; A(i, 1:21)];
    end
end


for i = 1:totalcase
    quadCoeff{i} = LHS{i}\RHS{i};
    quadCoeff{i}
end


for i = 1: m
    if isexplained(i) == 1
        continue;
    end
    
    temp = A(i, 1:5)';
    for j = 1:totalcase
        M = quadCoeff{j};
        if  max(abs(A(i, 6:21)' - flipba(M)' * temp)) <1e-3
            isexplained(i) = 1;
            break;
        end
    end
    if isexplained(i) == 0
        unexplained = [unexplained; A(i, 1:21)];

    else
        explained = [explained; A(i, 1:21)];

    end
end
explained = unique(explained, 'row');
unexplained = unique(unexplained, 'row');

fileid = fopen('temp.tex', 'w');

for caseno = 1:totalcase
    M = int64(quadCoeff{caseno});
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
    g
    gsym = str2sym(g);
    quadformulae(caseno) = gsym;
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
%subs(gsym, {k, x, y, z, w}, {5, -3, -1, -1, -1})


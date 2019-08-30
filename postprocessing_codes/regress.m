% inputsample
%b123 b234     b1234  b   b1  b2  b3  b4  ba  b1b2       b2b3        b1ba
% -8 -3 -3 -3  6:     0   0   0   8   3   0   0  -8  -3  -8  -3   0   8   8  -8  -3
% -8 -3 -3 -3  6:     0   0   8   0   3   0  -8   0  -3  -8   0  -3   8  -8   8  -3
%
%%

f = fopen("4varCasesNotDone.txt");

A = fscanf(f, "%f %f %f %f %f: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", [21, Inf]);
A = A';
totalcase = 16;
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
    
    % regress on some chosen cases
    
    if A(i, 6) == 0 && A(i, 7) == -x && A(i, 8) == 0 && A(i, 9) ==0 && A(i, 12)== x ...
            && A(i, 14) == z + w + k && A(i, 10) == -y  && A(i, 20) == -x - y + z && A(i, 16) == y
        caseno = 1;
        
    end
    
    if caseno ~= 0
        LHS{caseno} = [LHS{caseno}; A(i, 1:5)];
        RHS{caseno} = [RHS{caseno}; A(i, 6:21)];
        isexplained(i) = 1;
        %explained = [explained; A(i, 1:21)];
    end
end



quadCoeff{1} = LHS{1}\RHS{1}; 
quadCoeff{1}

X1 = matrixtoquad(quadCoeff{1});
for i = 2:totalcase
    quadformulae(i) = flipbits(2*(i-1), X1);
    quadformulae(i)
    quadCoeff{i} = quadtomatrix(quadformulae(i));
end


%%

batchsize = int64(20);
m = int64(m);
inittime=cputime;
for batchnum = 1:(idivide(m, batchsize)+1)
    parfor i = ((batchnum-1)*batchsize+1): min(m, batchnum*batchsize)
        temp = A(i, 1:5)';
        for j = 1:totalcase
            M = quadCoeff{j};
            if  max(abs(A(i, 6:21)' - flipba(M)' * temp)) <1e-3 ...
                    || max(abs(A(i, 6:21)' - M' * temp)) <1e-3
                isexplained(i) = j;
                break;
            end
        end
        if isexplained(i) == 0
            unexplained = [unexplained; A(i, 1:21)];
            fprintf("%d not done, ", i);
            
        else
            explained = [explained; A(i, 1:21)];
            fprintf("%d OK by %d, ", i, j);
        end
    end
    fprintf("\n %f %% done \n", double(batchnum)*double(batchsize)/double(m)*100)
end
explained = unique(explained, 'row');
unexplained = unique(unexplained, 'row');


%%

fileid = fopen('temp.tex', 'w');

for caseno = 1:totalcase
    syms b_1 b_2 b_3 b_4 b_a k x y z w
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

function X = flipba(Y)
X(:, 1) = Y(1:5, 1) + Y(1:5, 6);
X(:, 6) = -Y(:,6);
for j = 2:5 % b_(j-1) coeff
    X(:, j) = Y(:, j) + Y(:, 11+j);
    X(:, 11+j) = -Y(:, 11+j);
end
for j = 7:12
    X(:, j) = Y(:, j);
end

end
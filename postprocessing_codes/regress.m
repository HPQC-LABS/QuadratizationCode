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

    %manually find the correct grouping
    
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
    
    if caseno ~= 0 %ignore some quad because only need one quad per function
        LHS{caseno} = [LHS{caseno}; A(i, 1:5)];
        RHS{caseno} = [RHS{caseno}; A(i, 6:21)];
    end
end

%check whether have got all cases in the file
explained = zeros(0, 5);
for i = 1:8
    quadCoeff{i} = LHS{i}\RHS{i};
    quadCoeff{i}
    explained = [explained; LHS{i}];
end
explained = unique(explained, 'row');
allcase = unique(A(:, 1:5), 'row');

isequal(explained, allcase)

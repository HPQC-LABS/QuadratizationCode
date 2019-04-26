function genInputFolders(X, N, filehasno) % X the number cores, N the range of coeff, file has no. if fileno = 1
tic
if nargin == 0
    X = 500;
    N = 1;
    filehasno = 0;
end
allcases = zeros(0, 16);
a = zeros(1, 16);
for k = 0: int32((2*N+1)^16)
    if mod(k, 1000000) == 0
        toc
    end
    %a(1..10) are cubic coeff a123, a124, a125, a134, ..., a345; 
    %a(11..15) are quartic coeff a1234, a1235, a1245, a1345, a2345
    for j = int32(1) : 16
        a(j) = mod(k/(2*N+1)^(j-1), 2*N+1) - N;
    end
    
    %use symmetry 
    increasing = 1;
    for j = 11:14
        if a(j) > a(j+1)
            increasing = 0;
            break;
        end
    end
    var5 = 1;
    % check whether b1 appear 
    if a(1) == 0 && a(2) == 0 && a(3) == 0 && a(4) == 0 && a(5) == 0 && a(6) == 0 && a(11) == 0 && a(12) == 0 && a(13) == 0 && a(14) == 0
        var5 = 0;
    end
    % check whether b2 appear
    if a(1) == 0 && a(2) == 0 && a(3) == 0 && a(7) == 0 && a(8) == 0 && a(9) == 0 && a(11) == 0 && a(12) == 0 && a(13) == 0 && a(15) == 0
        var5 = 0;
    end
    
    if a(1) == 0 && a(4) == 0 && a(5) == 0 && a(7) == 0 && a(8) == 0 && a(10) == 0 && a(11) == 0 && a(12) == 0 && a(14) == 0 && a(15) == 0
        var5 = 0;
    end
    
    if a(2) == 0 && a(4) == 0 && a(6) == 0 && a(7) == 0 && a(9) == 0 && a(10) == 0 && a(11) == 0 && a(13) == 0 && a(14) == 0 && a(15) == 0
        var5 = 0;
    end
    
    if a(3) == 0 && a(5) == 0 && a(6) == 0 && a(8) == 0 && a(9) == 0 && a(10) == 0 && a(12) == 0 && a(13) == 0 && a(14) == 0 && a(15) == 0
        var5 = 0;
    end
    if increasing == 1 && var5 == 1
        allcases = [allcases; a];
    end
end
l = size(allcases, 1);

for fileno = 1 : X
    mkdir(sprintf('job_%d', fileno));
    cd(sprintf('job_%d', fileno));
    if filehasno == 1
        f = fopen(sprintf('input_coeff_%d.txt', fileno),'w');
    else
        f = fopen('input_coeff.txt', 'w');
    end
    for k = 1 + floor((fileno-1) * l/ X) : floor(fileno * l/ X)
        fprintf(f, "0,0,0,0,0,0,0,%d,0,0,0,%d,0,%d,%d,%d,0,0,0,%d,0,%d,%d,%d,0,%d,%d,%d,%d,%d,%d,%d\n", allcases(k, 1), allcases(k,2), allcases(k,4), allcases(k,7), allcases(k,11), allcases(k,3), allcases(k,5), allcases(k,8), allcases(k,12), allcases(k,6), allcases(k,9), allcases(k,13), allcases(k,10), allcases(k,14), allcases(k,15), allcases(k,16));
    end  
    cd ..
end
toc
end

X = 500; % number of cores
N = 3; 

allcases = zeros(0, 6);
for a1 = -N: N
    for a2 = a1 : N
        for a3 = a2 : N
            for a4 = a3 : N
                for a5 = a4 : N
                    for a6 = -N : N
                        allcases = [allcases; a1, a2, a3, a4, a5, a6];
                    end 
                end
            end
        end
    end
end
l = size(allcases, 1);

for fileno = 1 : X
    mkdir(sprintf('job_%d', fileno));
    cd(sprintf('job_%d', fileno));
    f = fopen(sprintf('input_coeff_%d.txt', fileno),'w');
    for k = 1 + floor((fileno-1) * l/ X) : floor(fileno * l/ X)
        for m = 1: 5
            fprintf(f, "%d,", allcases(k, m));
        end
        fprintf(f, "%d\n", allcases(k, 6));
    end  
    cd ..
end


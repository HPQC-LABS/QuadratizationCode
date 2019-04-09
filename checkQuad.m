function checkQuad(s)
% call checkQuad('has_quadratization.txt')

fileID = fopen(s, 'r');
temp = fgetl(fileID);

A = sscanf(temp, "%f");
a1234 = A(3);
a2345 = A(4);
a3451 = A(5);
a4512 = A(6);
a5123 = A(7);
a12345 = A(8);

n = 6;
allCombos=dec2bin(0:2^n-1) -'0';
b=mat2cell(allCombos,2^n,ones(1,n));
c=ones(2^n,1);
% columns are in order specified by Gruber
% columns are 1 b1 b2 ... b5 ba b1b2 b1b3 ... b4b5  b1ba b2ba ... b5ba
allbits=c;
for i = 1 : n
    allbits=[allbits b{i}];
end
for i=1:n-1
    for j=i+1:n-1
        allbits=[allbits b{i}.*b{j}];
    end
end
for i = 1: n-1
    allbits = [allbits b{i}.*b{n}];
end


while true
    temp = fgetl(fileID);
    if(temp == "$")
        break;
    end
    coeffsQ = sscanf(temp, "%f");
    
    LHS=a1234.*b{1}.*b{2}.*b{3}.*b{4} + a2345.*b{2}.*b{3}.*b{4}.*b{5} + a3451.*b{3}.*b{4}.*b{5}.*b{1} ...
        + a4512 .* b{4} .*b{5}.* b{1} .* b{2} + a5123 .* b{5}.* b{1} .* b{2} .* b{3} ...
        + a12345 .*b{1}.*b{2}.*b{3}.*b{4}.*b{5};
    RHS = allbits*coeffsQ;
    RHS = min(RHS(1:2:2^n -1), RHS(2:2:2^n));
    
    if max(abs(LHS(1:2:2^n-1)-RHS))<1e-5
        fprintf("yes this is a quadratisation\n")
        
    else
        fprintf("not a quadratisation\n");
    end
end

end

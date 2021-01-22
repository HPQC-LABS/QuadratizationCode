
n=4;

allCombos = zeros(0,4);
for i1 = [-1 0 1]
    for i2 = [-1 0 1]
        for i3 = [-1 0 1]
            for i4 = [-1 0 1]
                allCombos=[allCombos;i1 i2 i3 i4];
            end
        end
    end
end


t1 = allCombos(:,1);
t2 = allCombos(:,2);
t3 = allCombos(:,3);
ta = allCombos(:,4);
c = ones(3^n,1);
   
%%
allbits = [t1.*t2 t1.*t3 t1.*ta t2.*t3 t2.*ta t3.*ta t1 t2 t3 ta c];

warning('off', 'MATLAB:rankDeficientMatrix');

coeffsQ = ones(11,1);
numbers = int64(1:3:3^n-2);

N = 4;
quad = zeros(2*N+1,11);
%%

for alpha123 = [-N:-1 1:N]
    LHS = alpha123.*t1.*t2.*t3;
    LHS = LHS(1:3:3^n-2);
    
    
for i = int64(0): int64(2^(3^(n-1)-1))
    indexlist1 = bitget(i, 3^(n-1):-1:1);
    for j = int64(0):int64(2^(3^(n-1)-1))
        indexlist2 = bitget(j, 3^(n-1):-1:1);
        indexlist = indexlist1 + indexlist2 + numbers;
        A = allbits(indexlist,:);                     
        coeffsQ = A\LHS;           
        RHS = allbits*coeffsQ;                              
        RHS = min(reshape(RHS,3,[]))';


       if max(abs(LHS-RHS))<1e-5 
           quad(1+N+alpha123, :)=coeffsQ; %quadratization coefficients stored here
           fprintf("found quad for %dt1t2t3 = %+.1ft1t2 + %+.1ft1t3 + %+.1ft1ta + %+.1ft2t3 + %+.1ft2ta + %+.1ft3ta + %+.1ft1 + %+.1ft2 + %+.1ft3 + %+.1fta + %+.1f", alpha123, coeffs(:,1),coeffs(:,2),coeffs(:,3),coeffs(:,4),coeffs(:,5),coeffs(:,6),coeffs(:,7),coeffs(:,8),coeffs(:,9),coeffs(:,10),coeffs(:,11));
       end          
       if max(abs(quad(1+N+alpha123,:))) > 1e-5
             break;
       end
       
    end
end

end







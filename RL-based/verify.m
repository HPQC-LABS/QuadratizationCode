function [verified,const_terms] = verify(cache, coef)
	
    coeffs_size = cache.coeffs_size;
    
    coef = reshape(coef, [], coeffs_size);
    allbits = reshape(cache.allbits, 2^cache.n, coeffs_size);
    LHS = reshape(cache.LHS, 1, []);
    
    RHS = coef*allbits';
    if cache.aux
        RHS = min( RHS(:,1:2:end), RHS(:,2:2:end) ); % when using aux
    end

    const_terms = -min(RHS, [], 2);
    RHS = RHS + const_terms + min(LHS);
    
    % checking if every single state is equal to the min
    % of the corresponding states given by the COMIGs
    if min(RHS) == LHS
        verified = 1;
        fprintf('Hooray!!!\n');
    else
        verified = 0;
        fprintf('Try Again\n');
    end
end
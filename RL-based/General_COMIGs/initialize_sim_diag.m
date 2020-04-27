function [Q,IE] = initialize_sim_diag(LHS)
	Q = eye(size(LHS));
	IE = [ones(1,length(Q)-1),0];

	A = LHS;
	A = Q'*A*Q;
	epsilon = 10*max(size(A))*eps(normest(A));
	for j=find((IE==1).*[1,IE(1:end-1)==0])	% unaffected by 'IE' update
		k = j+find(IE(j:end)==0,1,'first')-1;
		[V,d] = schur(A(j:k,j:k)); d = diag(d);
		% handle case d = [1, i, (1-eps)*i]
		[~,IS] = sort(roundn10(d,round(log10(epsilon))+3));
		d = d(IS).';
		Q(:,j:k) = Q(:,j:k)*V(:,IS);
		IE(j:k) = [abs(d(1:end-1)-d(2:end)) < epsilon,0];
	end
end

function x = roundn10(x,n)

factors = 10^(fix(-n));
x = round(x*factors)/factors;
end

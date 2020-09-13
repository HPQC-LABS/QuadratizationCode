function param = initialize_sim_diag(varargin)
	if nargin == 1
		A = varargin{1};
		Q = eye(size(A));
		IE = [ones(1,length(Q)-1),0];
		J = find((IE==1).*[1,IE(1:end-1)==0]);
		K = [];
		for j = J, K = [K, j+find(IE(j:end)==0,1,'first')-1 ]; end
		epsilon = 10*max(size(A))*eps(normest(A));
		cutoff  = -round(log10(epsilon))-3;
	else
		Q  = varargin{1}.Q;
		IE = varargin{1}.IE;
		J  = varargin{1}.J;
		K  = varargin{1}.K;
		epsilon = varargin{1}.epsilon;
		cutoff  = varargin{1}.cutoff;
		A = varargin{2};
		A = Q'*A*Q;
	end
	
	for i = 1:numel(J)
		j = J(i); k = K(i);
		[V,d] = schur(A(j:k,j:k)); d = diag(d);
		[~,IS] = sort(round(d,cutoff)); % handle case d = [1, i, (1-eps)*i]
		d = d(IS).';
		Q(:,j:k) = Q(:,j:k)*V(:,IS);
		IE(j:k) = [abs(d(1:end-1)-d(2:end)) < epsilon,0];
	end
	param.Q  = Q;
	param.IE = IE;
	param.J  = find((IE==1).*[1,IE(1:end-1)==0]);
	param.K  = [];
	for j = param.J, param.K = [param.K, j+find(IE(j:end)==0,1,'first')-1 ]; end
	param.epsilon = epsilon;
	param.cutoff  = cutoff;
	if nargin == 1
		param.d = real(diag(Q'*A*Q));
	else
		param.d = real(diag(Q'*varargin{2}*Q));
	end
end
%{
function x = roundn10(x,n)

factors = 10^(fix(-n));
x = round(x*factors)/factors;
end
%}

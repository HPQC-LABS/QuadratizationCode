function flag = commutes(RHS)
	flag = all (cellfun(@(x) all( all( abs( RHS{end}*x - x*RHS{end} ) < 10e-5 ) ), RHS(1:end-1)) );
end
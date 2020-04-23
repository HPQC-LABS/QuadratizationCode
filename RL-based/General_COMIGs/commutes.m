function flag = commutes(RHS)
	flag = all(all( abs( RHS{1} * RHS{2} - RHS{2} * RHS{1} ) < 10e-5));
end

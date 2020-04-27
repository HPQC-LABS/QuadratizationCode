for i = -2:2
	if i == 0, continue, end
	for j = i:2
		if (j==-2) || ( i == -2 && j == -1 ) || (j==0), continue, end
		for k = -2:2
			quantum_envelopes_brute_force([i, j, k]);
		end
	end
end

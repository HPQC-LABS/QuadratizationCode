function plot_nullspace_sizes

	factorizable = true;
	
	plot;
	plot( factorizable);
	plot(~factorizable);
end

function plot(factorizable)

	figure;
	for i = 1:10
		if i<=5
			subplot(5,2,2*i-1);
		else
			subplot(5,2,2*i-10);
		end
		file = ['Polynomials/polynomialsDeg-3_Terms-', int2str(i),'.mat'];
		load(file, 'sizes');
		if nargin == 1
			load(file, 'factorizable_mask');
			if factorizable == true
				sizes = sizes(factorizable_mask);
			else
				sizes = sizes(~factorizable_mask);
			end
		end
		histogram(sizes, 'BinMethod', 'integers')
		xlim([1,19]);
		xticks(1:2:19);
		%set(gca, 'YScale', 'log');
		%ylim([1,10e6]);
		title([int2str(i), ' Terms']);
		legend(sprintf('mean = %.1f, std = %.1f', mean(sizes), std(double(sizes)) ));
		if i == 1,	title('1 Term');	end
	end
end

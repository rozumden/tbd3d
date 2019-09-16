function tuple = nthtuple(n, k)
% Intended for relatively fast egenration of unique k-tuples - n-th tuple for each 'n' when tuples are 'ordered'. Usage - generate 'n' randomly and uniquely, then convert 'n' to tuples without having to check for uniqueness etc
% The method is nor very nice and maybe not very efficient - I couldn't come up with a better one though. It relies on correct rounding of doubles to integers (after square roots etc) so it may break down for some (large) combinations of 'n' and 'k'
% Tested and works for
%	1/ k=2, n up to 1e9 (corresponds to all pairs in > 40k points)
% 	2/ k=4, n up to >1e9 (corresponds to all quadruples in 400 points)

n = n(:);
tuple = zeros(numel(n), k);

for k=k:-1:1
	switch(k)
		case 1
			tuple = n;
		case 2
			% solve directly in quadratic case
			tuple(:,2) = ceil((1+sqrt(1+8*n))/2);
			tuple(:,1) = n-(tuple(:,2)-1).*(tuple(:,2)-2)/2; % shortcut to k=1 case
			return;
		otherwise
			% solve for last entry via approximate solution to k-th order eq and do the rest recursively
			tuple(:,k) = ceil(nthroot(n*factorial(k), k)+(k-1)/2); % lower bound for the correct solution (the correct solution can be this+1 (tested for k=4))
			val1 = tuple(:,k)-1;
			val2 = tuple(:,k);
			for i=1:k-1
				val1 = val1.*(tuple(:,k)-1-i)/(i+1); % computes nchoosek(tuple(k)-1,k) (defined as 0 for n<k)
				val2 = val2.*(tuple(:,k)-i)/(i+1); % computes nchoosek(tuple(k),k)
			end 
			where = val2 < n; tuple(where,k) = tuple(where,k) + 1; val1(where) = val2(where); % note: for same large 'k' and large 'n' this may fail and may actually need a loop adding +1 (maybe not - this requires further testing for large input)
			%fprintf('k=%d, #where=%d/%d\n', k, nnz(where), numel(where)); % DBG
			% recursive call
			n = n - val1;
			% tuple(:,1:k-1) = nthtuple(n-val1, k-1);
	end
end
end

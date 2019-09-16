function coeffs = parabola4p(pts)
% gives implicit equation of the parabola that passes through given 4 points ("minimal problem"), in 2D. In general 2 parabolas pass though 4 pts, so these two solutions are given.
% Does not return degenerate solutions (complex, straight lines etc)
% 
% pts - 4x2 array of point coords
% coeffs - 2x5 array of coeffs for 2 solutions. In each row: [c1 c2 c3 c4 c5] in the parabola equation (c1*x+c2*y)^2+c3*x+c4*y +c5 = 0 where [x y] are first/second resp coordinate (column) in pts.
% 
% note: see https://www.mathpages.com/home/kmath546/kmath546.htm for notation etc

eps = 1e-10; % zero

% transform coords so that p4=0 (simplifies eqs, offset will be added to the result)
offset = pts(4,:);
pts = pts(1:3,:) - offset;

% common relations
v = cross(pts(:,1), pts(:,2)).';
d = -prod(v)*sum(v); % discriminant of the quadratic eq for A
if(d < eps) coeffs = []; return; end % obviously degenerate solution - three pts in a line or 1 pt inside triangle
num = v*prod(pts,2); % numerator term for the A solution (equals 'b'/2 term of the quad eq)
[~, ix] = max(abs(v)); % index of the max 'v' - for best numerical stability in the division later

% try parabola as (Ax+y)^2
[coeffs count] = xproblem(pts, offset);

if(count < 2) % single solution (unlucky orientation) - flip xy and continue
	v = -v; num = -num; % flipping xy changes polarity of 'v'
	[coeffs2 count] = xproblem(pts(:,[2 1]), offset([2 1])); % parabola as (x+Cy)^2
	if(count < 2) % super unlucky situation, two parabolas with perpendicular axes, parallel to coordinate system -> combine
		coeffs = [coeffs; coeffs2(:, [2 1 4 3 5])]; % flip xy of the second solution
	else % keep second solution only, easy way how to get rid of duplicates
		coeffs = coeffs2(:, [2 1 4 3 5]); % flip xy
	end
end

function [coeffs count] = xproblem(pts, offset)
% finds coeffs of the parabola (Ax+y)^2+Dx+Ey+F = 0, ie parabola passes through the (given three pts and origin)+offset and y^2 coeff is fixed to 1 (cannot be y=x^2 etc)
% pts - 3x2 for three points as [x(:) y(:)] coords; offset 1 point
%
% returns coeffs as 2x5 array (or less) of coeffs [A 1 D E F] for parabola (Ax+y)^2+Dx+Ey+F=0
% count - how many theoretical solutions are found (some can be degenerate and removed)

% quad coeff solutions (A)
denom = v*pts(:,1).^2; % denominator of the solution (equals 'a' term of the quad eq)
if(abs(denom) < eps) % linear eq, single solution
	A = -(v*pts(:,2).^2)/(2*num); % note: num = 0 with positive discriminant cannot happen
	count = 1;
else % quadratic eq, two solutions
	A = (-num + [-1 1]*sqrt(d))/denom; % x coeff solutions
	count = 2;
end

% linear coeffs (D,E)
u2 = (A.*pts(:,1) + pts(:,2)).^2; % (A*x+y)
temp = cross(repmat(pts(:,2), [1 size(u2,2)]), u2, 1); D = temp(ix, :)/v(ix);
temp = cross(u2, repmat(pts(:,1), [1 size(u2,2)]), 1); E = temp(ix, :)/v(ix);

% remove degenerate solution (quadratic eq can describe two parallel line, as in (x-1)^2=1 etc)
% the condition is D!=AE, follows from expanding 'distance to straight line' (ax+by+x)^2=const
keep = abs(D-A.*E) > eps;
A = A(keep); D = D(keep); E = E(keep);

% apply offset
coeffs = [A(:) ones(numel(A),1) D(:)-2*(offset(1)*A(:).^2+offset(2)*A(:)) E(:)-2*(offset(1)*A(:)+offset(2)) A(:).^2*offset(1)^2+offset(2)^2+2*A(:)*offset(1)*offset(2)-D(:)*offset(1)-E(:)*offset(2)]; % a little bit complicated set of coeffs which follows from writing the parabola in matrix form as xT*m*x+bT*x+c=0 and subs (x-offset)
end
end
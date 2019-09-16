function [m f] = project2pyramid(m, f, eps)
% Subproblem in estimateFM_* - projection of (m,f) values to feasible "pyramid"-like intersection of convex sets (roughly all positive, m<=1, m>=f)
% Fast version using Dykstra's algorithm, see eg Ryan J. Tibshirani, Dykstra’s Algorithm, ADMM, and Coordinate Descent: Connections, Insights, and Extensions, NIPS 2017
% Note: the sets are 1/ f*>0, m>0, m<1 and 2,3,4,../ m - f_i*eps > 0 for i=RGB... (ie 4 conv sets for RGB image f)
%
% m - single column vector
% f - RGB columns vector-triplet (vec3) (or arbitrary dimension shaped likewise)
% eps - see estimateFM_3map_*; inverse slope of the f<=m/eps constraint for each channel. eps=0 means no constraint (only m in [0,1], f>0), eps=1 means f<=m etc.
% NOTE: using this function for eps=0 makes no sense, such projection can be done in one step in the calling function

maxiter = 10; % number of whole cycles, 0=only positivity

u = [m f];
d = size(u,2); % number of conv sets
z = cell(1,d);
[z{:}] = deal(zeros(size(u))); % auxiliary vars (sth like projection residuals)
n = [-1 eps]./sqrt(1+eps^2); % dividing plane normal vector (for projection to oblique planes)

for iter=1:d*maxiter+1 % always end with projection to the first set
	u_old = u;

	idx = mod(iter-1,d)+1; % set index

	% projection
	if(idx == 1) % projection to f>0, 0<m<1
		u = u + z{1};
		u(u<0) = 0; % f,m > 0
		u(u(:,1) > 1,1) = 1; % m < 1
	else % one of the oblique sets
		% projection
		u = u + z{idx};
		w = u(:,idx)*eps > u(:,1); % points outside of C_idx
		proj = u(w,[1 idx])*n.'; % projection to normal direction
		u(w,[1 idx]) = u(w,[1 idx]) - n.*proj;
	end

	% auxiliaries
	z{idx} = u_old + z{idx} - u;
end

m = u(:,1);
f = u(:,2:end);
end
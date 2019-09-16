function res = IF(condition, true_action, false_action)
% substitutes C++ like ? operator. When either of actions can be evaluated only if condition holds/doesn't hold, it needs to be passed as function handle.
%
% Usage - simple: x = IF(a>b, 1, 0);
% Usage - deferred evaluation: x = IF(~isempty(a), @()a(1), 0); % a(1) is evaluated only if `a` is not empty

if(condition)
	res = true_action();
else
	res = false_action();
end
end

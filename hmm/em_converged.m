function [converged, decrease] = ...
	em_converged(ll1, ll0, threshold, check_increased)

if nargin < 3, threshold = 1e-10; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
  if abs(ll1 - ll0) > 0.5*abs(ll1) && ~isinf(ll0)  % allow for a little imprecision
    %fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', ll0, ll1);
    decrease = 1;
	converged = 1;
	return;
  end
end

delta_loglik = abs(ll1 - ll0);
avg_loglik = (abs(ll1) + abs(ll0) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end

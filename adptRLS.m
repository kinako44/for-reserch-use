[w, e] = adptRLS(x, desired, tap, lambda, delta)
% RLSによる線形適応フィルタ
% x			: input signal
% desired	: desired signal
% tap 		: length of filter
% lambda	: RLS forgetting factor
% delta		: small positive constant
%
% w 		: filter coefficients
% e 		: error

% Parameters
iter = min(length(x), length(desired));

% Initialize
e = zeros(iter, 1);
w = zeros(tap, 1);
P = delta * eye(tap);
w = zeros(tap, 1);
xi = zeros(tap, 1);
	
% Execution
for j = 1 : iter
	xi = [x(j) ; xi(1:end-1)];
	g = (P * xi) / (lambda + xi.' * P * xi);
	v = desired(j) - w.' * xi;
	w = w + g * v;
	P = 1/lambda * (P - g* xi.' * P);
	e(j) =  v;
end
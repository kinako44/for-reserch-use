function [w, error] = adptNLMS(x, desired, tap, mu0)
% NLMSアルゴリズムに基づいた線形適応フィルタ

% x			: 入力信号 
% desired	: 所望信号
% tap		: 適応フィルタのタップ数
% mu0		: ステップサイズ
%
% w			: 適応フィルタの係数ベクトル
% error		: 誤差ベクトル


% Parameters
iter = length(x);		% iterations

% Initialize
error = zeros(iter, 1);
w = zeros(tap, 1);
xi = zeros(tap, 1);

% Computing
for j = 1 : iter
	xi = [x(j) ; xi(1:end-1)];
	y = xi.' * w;
	error(j) = desired(j) - y;
	mu = mu0 / (xi.' * xi);		% normalize
	w = w + mu * error(j) * xi;		% update coefficients
end

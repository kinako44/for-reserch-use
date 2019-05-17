function [ y ] = calcFLANN( x, weight )
% FLANNの出力信号を計算する関数

% x 		: 入力信号
% weight	: 重み（adptFLANNで同定したものを想定しています）
%			  (weightはN*(2P+1)の行列)
%
% y 		: 出力信号


% parameter
P = (size(weight, 2) - 1) / 2;

% execution
y = filter(weight(:, 1), 1, xt);
for i = 1 : P
	y = y + filter(weight(:, 2*i), 1, sin(i*pi*x)) + filter(weight(:, 2*i+1), 1, cos(i*pi*x));
end


end
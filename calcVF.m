function [ y ] = calcVF( x, kernel )
% Volterraフィルタの出力を計算
%
% x 		: 入力信号
% kernel	: Volterra核(adptVF2, adptVF3に対応．1*q cell array)
%
% y			: 出力信号


% parameter
iter = length(x);
q = length(kernel);
tap = length(kernel{1});

% execution
y = filter(kernel{1}, 1, x);
if (q > 1)
	xvec = zeros(tap, 1);
	for i = 1:iter
		xvec = [x(i) ; xvec(1:end-1)];
		y(i) = y(i) + xvec.' * kernel{2} * xvec;
		if (q > 2)
			for j = 1:tap
				y(i) = y(i) + xvec(j) * xvec.' * kernel{3}(:,:,j) * xvec;
			end
		end
	end
end
end

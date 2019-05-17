function [ weight, error ] = adptFLANN( x, desired, memory, P, mu0 )
% FLANNを用いた非線形システム同定用のfunction
% 更新アルゴリズムはLMS, 多項式展開は三角多項式
% 繰り返し回数はxおよびdesiredの長さと同じです

% <Input>
% x			: 入力信号 (N*1 vector)
% desired	: 所望信号 (N*1 vector)
% memory	: 記憶長 (scalar)
% P			: 関数展開次数 (scalar)
% mu		: ステップサイズ (scalar)
%
% <Output>
% weight	: 重み (memory*2P+1 matrix)
% error		: 誤差信号 (N*1 vector)

if (memory < 0 || P < 0 || mu0 < 0)
	fprintf('error\n');
	weight = 0;
	error = 0;
	return;
end
if (length(x)~=length(desired))
	fprintf('ベクトルの長さが異なります\n');
	return;
end

%% Parameters
iter = length(x);
weight = zeros(memory, 2*P+1);
error = zeros(iter, 1);
expanded = zeros(iter, 2*P);
exvec = zeros(memory, 2*P+1);


% calculate functional expansion signal
for j = 1:P
	expanded(:, 2*j-1:2*j) = [sin(j*pi*x) cos(j*pi*x)];
end

%% Exection
for i = 1:iter
	exvec = [x(i) expanded(i, :) ; exvec(1:end-1, :)];
	
	% Output signal
	y = 0;
	for j = 1 : 2*P+1
		y = y + weight(:, j).' * exvec(:, j);
	end

	error(i) = desired(i) - y;
	
	% Update weight
	for j = 1 : 2*P+1
		weight(:, j) = weight(:, j) + mu0 * error(i) * exvec(:, j);
	end
	
	
end

end

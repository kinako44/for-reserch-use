function [ kernel, error ] = adptVF3( x, desired, tap, q, mu )
% LMSによる適応Volterraフィルタ
% VF3は2次以降の所望信号に誤差信号を用います(コメント参照)
% ステップサイズとしてmuを1つのみ与えていますが，次数によって値を変えたほうがいいかと思います
%
%
% <Input>
% x			: input signal (N*1 vector)
% desired	: desired signal (N*1 vector)
% tap		: memory size of volterra filter (scalar)
% q			: order of volterra filter (scalar)
% mu		: step size parameter (scalar)
%
% <Output>
% kernel	: identified volterra kernel(1*q cell)
% error		: error signal (N*q matrix)


%% Parameters
iter = length(x);

kernel = cell(1, q);
% fix later
kernel{1} = zeros(tap, 1);
if (q > 1)
	kernel{2} = zeros(tap, tap);
	if (q > 2)
		kernel{3} = zeros(tap, tap, tap);
		if (q > 3)
			disp('error : Order > 3');
			kernel = 0;
			error = 0;
			return;
		end
	end
end

%% Execution
% init
error = zeros(iter, q);
xvec = zeros(tap, 1);
mu1 = mu;
mu2 = 0.5*mu;		% 調整用
mu3 = 0.5*mu;		% 調整用

for i = 1:iter
	% Linear
	xvec = [x(i) ; xvec(1:end-1)];
	y = xvec.' * kernel{1};
	error(i, 1) = desired(i) - y;
	kernel{1} = kernel{1} + mu1 * error(i, 1) * xvec;
end

% 2nd order	
if (q > 1)
	xvec = zeros(tap, 1);
	k = zeros(tap, tap);
	for i = 1:iter
		xvec = [x(i) ; xvec(1:end-1)];
		y = xvec.' * k * xvec;
		error(i, 2) = error(i, 1) - y;			% 1次Volterra核の計算で得られた誤差信号を所望信号として利用
		k = k + mu2 * error(i, 2) * (xvec * xvec.');
% 		k = k + mu2 * error(i, 2) * (xvec * xvec.');
	end
	kernel{2} = k;
end

% 3rd order
if (q > 2)
	xvec = zeros(tap, 1);
	k3 = zeros(tap, tap, tap);
	for i = 1:iter
		xvec = [x(i) ; xvec(1:end-1)];
		y = 0;
		xmat = xvec * xvec.';
		for j = 1:tap
			y = y + xvec(j) * sum(sum(k3(:,:,j) .* xmat));
		end
		error(i, 3) = error(i, 2) - y;			% 2次Volterra核の計算で得られた誤差信号を所望信号として利用
		
		for j = 1:tap
			k3(:,:,j) = k3(:,:,j) + mu * error(i, 3) * xvec(j) * xmat;
		end
	end
	kernel{3} = k3;
end


end

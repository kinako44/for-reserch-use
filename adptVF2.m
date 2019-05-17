function [ kernel, error ] = adptVF2( x, desired, tap, q, mu )
% LMSによる適応Volterraフィルタ
% VF2はkernelの更新にすべて共通の誤差信号を用います
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


%% parameters
iter = length(x);
kernel = cell(1, q);

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
mu1 = mu*40;
mu2 = mu;

fprintf('processing 1st order\n');
k1 = zeros(tap, 1);
for i = 1:iter
	% Linear
	xvec = [x(i) ; xvec(1:end-1)];
	y = xvec.' * k1;
	error(i, 1) = desired(i) - y;
	k1 = k1 + mu1 * error(i, 1) * xvec;
	if (i > iter-1000)
		kernel{1} = kernel{1} + k1/1000;
	end
end


% 2nd order	
if (q > 1)
	fprintf('processing 2nd order\n');
	xvec = zeros(tap, 1);
	k2 = zeros(tap, tap);
	for i = 1:iter
		xvec = [x(i) ; xvec(1:end-1)];
		y = xvec.' * k2 * xvec;
		error(i, 2) = error(i, 1) - y;
		k2 = k2 + mu2 * error(i, 2) * (xvec * xvec.');
		if (i > iter-1000)
			kernel{2} = kernel{2} + k2 / 1000;
		end
	end
% 	kernel{2} = k2;
end

% 3rd order
if (q > 2)
	fprintf('processing 3rd order\n');
	xvec = zeros(tap, 1);
	k3 = zeros(tap, tap, tap);
	for i = 1:iter
		xvec = [x(i) ; xvec(1:end-1)];
		y = 0;
		xmat = xvec * xvec.';
		for j = 1:tap
			y = y + xvec(j) * sum(sum(k3(:,:,j) .* xmat));
% 			y = y + xvec(j) * xvec.' * kernel{3}(:,:,j) * xvec;
		end
		error(i, 3) = error(i, 2) - y;
		for j = 1:tap
			k3(:,:,j) = k3(:,:,j) + mu * error(i, 3) * xvec(j) * xmat;
		end
	end
	kernel{3} = K3;
end

error = error(:, 2);
end

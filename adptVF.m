function [ kernel, error ] = adptVF( x, desired, tap, q, mu )
% LMSによる適応Volterraフィルタ
% VFはkernelの更新に共通の誤差信号を用います
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
% error		: error signal (N*1 vector)


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
error = zeros(iter, 1);
xvec = zeros(tap, 1);
mu1 = mu;
mu2 = mu*0.5;	% 調整用
mu3 = mu*0.5;


for i = 1:iter
	% output signal
	xvec = [x(i) ; xvec(1:end-1)];
	y = xvec.' * kernel{1};
	if (q > 1)
		y = y + xvec.' * kernel{2} * xvec;
		if (q > 2)
			for j = 1:tap
				y = y + xvec(j) * xvec.' * kernel{3}(:,:,j) * xvec;
			end
		end
	end
	
	error(i) = desired(i) - y;

	% update kernel
	kernel{1} = kernel{1} + mu1 * error(i, 1) * xvec;
	if (q > 1)
		kernel{2} = kernel{2} + mu2 * error(i, 2) * (xvec * xvec.');
		if (q > 2)	
			for j = 1:tap
				kernel{3}(:,:,j) = kernel{3}(:,:,j) + mu * error(i, 3) * xvec(j)*(xvec * xvec.');
			end
		end
	end

end

end

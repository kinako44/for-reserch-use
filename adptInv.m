function h_inv = adptInv(h, iter, N, mu)
% 線形逆フィルタ設計スクリプト
% based on Filtered-x LMS

% h : 入力インパルス応答
% iter : 繰り返し回数
% N : 逆フィルタのタップ数
% mu : stepサイズ
%
% h_inv : 逆フィルタ
%

% 逆フィルタのタップ数は，入力したインパルス応答の4倍程度あれば
% だいたい同定できるかと思います．
% うまくいかない場合は，インパルス応答を正規化して実行してみてください．
% その場合は，設計した後で逆フィルタの係数を調整する必要があります．
%


if(size(h, 1) == 1)
	h = h.';
end

%% Setting
Nt = 1;							% 試行回数（平均をとる） 
Nw = 0.5;						% 学習信号の分散

% parameter
mu0 = 0.001;
delay = N / 2;							% sample

% results
me = zeros(iter, 1);			% 誤差
mh = zeros(N, 1);				% 逆フィルタの係数ベクトル


%% Execution
for j = 1 : Nt
	
	%% 学習信号の生成
	rng('shuffle');							% ランダマイザのリセット（ver.によって実行不可の場合コメントアウト）
	x = Nw * randn(iter, 1);
	%% Initialize
	h_inv = zeros(N, 1);
	ref = zeros(N, 1);
	input = zeros(N + length(h) - 1, 1);
	desired = zeros(iter, 1);
	desired(delay+1:end) = x(1:end-delay);
	e = zeros(iter, 1);
	
	%% Identification
	for i = 1 : iter
		input = [x(i) ; input(1:end-1)];
		
		% Calculate filtered reference signal
		ref = [h.' * input(1:length(h)) ; ref(1:end-1)];
		
		% Output signal
		tf = conv(h, h_inv);
		y = tf.' * input;
		
		error = desired(i) - y;
		
		% update coefficents
% 		mu = mu0 / (input.' * input);		% コメント解除でNLMSによる更新に変更
		h_inv = h_inv + mu * error * ref;
		
		e(i) = error;
	end
	%% Calculate means
	me = me + e / Nt;
	mh = mh + h_inv / Nt; 
end

% MSEの表示
figure,
plot(10 * log10(me.^2)),
title('MSE of design of the inverse filter'),
xlabel('Number of iterations'), ylabel('MSE (dB)');

% hvtoolによる特性の表示
hfvt = fvtool(h, 1, h_inv, 1, conv(h, h_inv));
legend(hfvt, 'H', 'H^{-1}', 'H\astH^{-1}');

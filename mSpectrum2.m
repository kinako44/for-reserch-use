function mSpectrum2( data, fs, varargin )
% 振幅スペクトルを表示するfunction
% 振幅1の信号を0dBで表示

% data			: N*1 array
%				: length(data) = 2^nにすることをおすすめします
%
% fs			: sampling frequency (Hz)
% varargin(1) 	: option (窓関数の種類を選択可能．デフォルトはblackman窓)
%				: 0 - ハン窓，1 - ハミング窓，2 - ブラックマン窓，3 - 矩形窓


% サイズ確認用
if (size(data, 1) < size(data, 2))
	data = data.';
end

N = length(data);		% No. of point

% optionの確認
if (isempty(varargin))
	window_type = 3;
else
	window_type = varargin{1, 1};
end

% 窓関数の選択
switch window_type
	case 0
		window = hann(N);
	case 1
		window = hamming(N);
	case 2
		window = blackman(N);	
	case 3
		window = ones(N, 1);
end

% 以下計算
xaxis = (0:N/2).' * fs / N;			% 範囲は0からfs
Y = fft(data.*window);
yaxis2 = abs(Y/N);					% パワーがpoint数に依存しないように
yaxis = yaxis2(1:N/2+1);			
yaxis(2:end-1) = 2*yaxis(2:end-1);	% 折り返し分のパワーを考慮

yaxis = 20*log10(yaxis);


%% 振幅特性
figure,
plot(xaxis, yaxis)
title('Amplitude characteristics')
xlabel('Frequency (Hz)', 'FontSize', 12)
ylabel('Amplitude (dB)', 'FontSize', 12)

%% 位相特性
% figure
% plot(xaxis, unwrap(angle(yaxis(1:N/2))))
% title('Phase characteristics')
% xlabel('Frequency (Hz)')
% ylabel('Phase (rad)');
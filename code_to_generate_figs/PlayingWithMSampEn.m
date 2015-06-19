N=3*10^4;


w = randn(N,1);
hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',N);
f = step(hcn);
Fs = 1;

% f1 = f/var(f);

% sampen(w, 2, 0.15);

x = w;


for tau=1:20
    y = 0;
    for j = 1:N/tau
        y(j) = (1/tau) * sum(x((j-1)*tau+1:j*tau));
    end
    e = sampen(y, 2, 0.15, 0, 0, 0);
    fprintf('tau=%d done.\n', tau);
    ySampEn(tau) = e(2);
end


%%
plot(white, 'xk');
hold on;
plot(pink, 'xm');
legend('White Noise', '1/f (Pink) Noise');
xlabel('Scale Factor, $\tau$');
ylabel('SampEn ($\tau$)');
CJWPlotV4('MSampEnWhitePink');

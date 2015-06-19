function [ ibi_interp ] = cjw_resp( x, fs )
%CJW_RESP Summary of this function goes here
%   Detailed explanation goes here

T = length(x)/fs;

br_pm_max = 45;
br_ps_max = br_pm_max/60;

br_pm_min = 16;
br_ps_min = br_pm_min/60;

dfs = 5;

dx = downsample(x, fs/dfs);

[b,a]=butter(4,[br_ps_min,br_ps_max]/(dfs/2),'bandpass');
y=filtfilt(b,a,dx);

t = 0:1/dfs:(length(y)-1)/dfs;

% plot(t, y);

% c is the vector of spikes at the crossings
c = diff(sign(y));

% r is the vector of spikes of magnitude 1
for i=1:length(c)
    if c(i)~=0
        r(i) = 1;
    else
        r(i) = 0;
    end
end


% br_ps_max is breaths per second (Hz)
% 5/br_ps_max is number of samples between breaths

% f is the vector of spikes with the minimum distance between each (based o
% n max breath rate)
last = -inf;
for i=1:length(r)
    if r(i)==1 && i-last > dfs/br_ps_max
        last = i;
        f(i) = 1;
    else
        f(i) = 0;
    end
end

% plot(c);
% hold on;
% scatter(find(f), zeros(size(find(f))));
% hold off

% times is the vector of times of the breaths
times = find(f);

ibi = diff(times);

ibi_fs = 4;

TQ=[1:1:T*ibi_fs]/ibi_fs;

ibi_interp = interp1(times(2:end)/dfs,ibi,TQ,'PCHIP');

end


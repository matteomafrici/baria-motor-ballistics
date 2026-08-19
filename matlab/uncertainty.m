function [a, Inc_a, n, Inc_n, R2] = uncertainty(p, rb)
%UNCERTAINTY  Ordinary least squares fit of Vieille's law in log-log space.
%   [a, Inc_a, n, Inc_n, R2] = uncertainty(p, rb)
%
%   Fits rb = a * p^n to the (p, rb) couples by linearising
%   ln(rb) = ln(a) + n * ln(p) and applying ordinary least squares.
%   Inc_a and Inc_n are the symmetric uncertainties on a and n.
%
%   Original file: Incertezze.m V 1.02, 04 August 2011 - Dossi

X = log(p);                          % log pressure
N = length(p);
Y = log(rb);                         % log burning rate

delta = N * sum(X.^2) - (sum(X))^2;  % denominator for q and m
m = (N * sum(X.*Y) - sum(X)*sum(Y)) / delta;   % slope = n
n = m;
q = (sum(X.^2)*sum(Y) - sum(X)*sum(X.*Y)) / delta;  % intercept = ln(a)

Y_eval = m .* X + q;                 % regression values
sig = (sum((Y - Y_eval).^2) / (N - 2))^0.5;  % RMS residual

Inc_n = sig * (N / delta)^0.5;       % uncertainty on n
Inc_q = sig * (sum(X.^2) / delta)^0.5;      % uncertainty on ln(a)

a = exp(q);                          % Vieille coefficient
a_new_up = a * exp(Inc_q);           % upper bound
a_new_down = a / exp(Inc_q);         % lower bound
Inc_a = (a_new_up - a_new_down) / 2; % symmetric uncertainty on a

M_Y = mean(Y);                       % regression R^2
dev_reg = sum((Y_eval - M_Y).^2);
dev_tot = sum((Y - M_Y).^2);
R2 = dev_reg / dev_tot;
end

function [Y,F]=spectralDecomp(y,plot)
% PURPOSE: Do a variance decomposition of data using spectral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
% y = vector (T x 1) with stationary time series
%
% plot = logical. If true, result is plotted. Default=false;
%
% Output:
%
% Y = vector (T/2+1 x 1) with periodogram of y
%
% F = vector (T/2+1 x 1) with frequencies on which Y is evaluated
%
% Usage:
%
% [Y,F]=spectralDecomp(y,varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Anders Thorsrud
% leifath@gmail.com
% 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


y=y(:);
T=numel(y);

Fs=T;
NFFT=T;

yf=fft(y,NFFT)/T;

Y=2*abs(yf(1:NFFT/2+1));
F=Fs/2*linspace(0,1,NFFT/2+1)./T;

if nargin==2
    % Plot single-sided amplitude spectrum.
    f=figure('name','Periodogram');
    h=plot(F,Y); 
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency')
    ylabel('|Y(f)|')
end;



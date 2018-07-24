% =======================================================================
% Determine lag order of AR(4) model with constant on simulated data using
% information criteria
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================
clearvars; clc; close all;
ENDO = xlsread('../data/AR4.xlsx',1);         % load data
pmax = 8;                                     % set maximum number of lags
const = 1;                                    % model with constant
nlag = LagOrderSelectionARp(ENDO,const,pmax); % compute criteria
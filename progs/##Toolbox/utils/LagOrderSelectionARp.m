clear all; close all; clc;
ENDO = xlsread('../data/data.xlsx','AR4');
% True AR order is 4
opt.maxlags = 10;
opt.const = 1;
nlag = LagOrderSelection(ENDO,opt);
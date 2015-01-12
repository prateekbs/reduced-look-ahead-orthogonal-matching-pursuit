function [ x, positions ] = generateSignal(TYPE, N, K  )
% Generate Sparse Signal/Source - Sooraj
%
% Inputs
%   TYPE	: type of the signal required
%               'StdGaussian' - Std Gaussian signals
%               'Bernoulli'  - all ones in  non-zero location
%               'BernoulliRndSign' - +1/-1 with equal probability
%               (Rademacher distribution)
%   N           : signal length
%   K           : sparsity level
% Outputs
%   x   	: generated signal
%   positions   : support set (locations in which signal is not zero)
%
% Purpose   	: To generate different sparse signals for simulations
% Reference :
%
% Author       	: Sooraj K. Ambat
% Email        	: sooraj@ece.iisc.ernet.in/sooraj.k.ambat@gmail.com
% Address      	: Ph.D. Scholar
%                 Statistical Signal Processing Lab
% 		  Department of Electrical Communication Engineering
%                 Indian Institute of Science, Bangalore, India-560 012.
%
% Date         	: 20 Aug 2011
% Last Modified  : 09 May 2012


% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Short Disclaimer: this script is for educational purpose only.

positions  = randperm(N);
positions = positions(1:K); % Location of non-zero rows in X
x = zeros(N,1);
if(strcmp(TYPE, 'StdGaussian'))
    x(positions,:) = randn(K,1);    % The k-row-sparse solution
elseif (strcmp(TYPE, 'BernoulliOne'))
    x(positions,:) = 1; 
elseif(strcmp(TYPE, 'BernoulliRndSign'))
    tmp = randi([0 1], K,1);    % The k-row-sparse solution (bernoulli distributed)    % The K-row-sparse solution
    tmp(tmp == 0) = -1;      %replace zeros with -1
    x(positions,:) = tmp; %randi([0 1], K,1);    % The k-row-sparse solution (bernoulli distributed)    % The k-row-sparse solution
else
    error('Unknown Option Specified to generate Signal');
end


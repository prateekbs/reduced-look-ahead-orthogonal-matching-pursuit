function [estX, estSupportSet] = SPwithInitialSupport(A, b, K1, initSupportSet)
% Subspace Pursuit algorithm implementation which supports partial 
% information about support-set
%

% Inputs
%   A           	: Sensing Matrix (size MxN)
%   b           	: Measurement Vector (size Mx1) b = Ax
%   K1           	: sparsity level
%   initSupportSet 	: initial Support-set
% Outputs
%   estX                : estimated sparse signal 'x'
%   estSupportSet     	: estimated support of x
%
% Purpose        	: To implement the Subspace Pursuit Algorithm with 
%			  initial Support.
% Reference     	: "Subspace Pursuit for Compressive Sensing Signal
%                          Reconstruction"-Wei Dai and Milenkovic, O.
%                          IEEE Transactions on Information Theory, 2009.
%                       
% Author            	: Sooraj K. Ambat
% Email              	: sooraj@ece.iisc.ernet.in/sooraj.k.ambat@gmail.com
% Address           	: Ph.D. Scholar,
%                         Statistical Signal Processing Lab,
%                         Indian Institute of Science, Bangalore, India-560 012.
% Date               	: 29 Sep 2011
% Modified       		: 27 May 2012

% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Short Disclaimer: this script is for educational purpose only.
%update residue with the initSupport
rk = b-A(:,initSupportSet)*(A(:,initSupportSet)\b);
currResNorm = norm(rk,2);
estSupportSet = initSupportSet;
[M,N] = size(A);
maxIter = M; %maximum number of iterations
for counter=1:maxIter
    prevestSupportSet = estSupportSet;    % storing prevestSupportSet value
    %matched filter choose best K1 support atoms
    matchFilter = abs(A'*rk);
    [~, index1] = sort(matchFilter, 'descend');
    interimSupportSet = union(estSupportSet, index1(1:K1)'); %K1<=|estSupportSet|<=2K1
    %choose best K1 using the interim estimate of signal
    interimEstX = zeros(N,1);
    interimEstX(interimSupportSet) = A(:,interimSupportSet)\b;
    [~,index2] = sort(abs(interimEstX), 'descend');
    estSupportSet = index2(1:K1)';
    %update residue
    rk = b-A(:,estSupportSet)*(A(:,estSupportSet)\b);
    prevResNorm = currResNorm;
    currResNorm = norm(rk,2);
    %check for halting criteria
    if(prevResNorm <= currResNorm)
        estSupportSet = prevestSupportSet;        
        break;
    end    
end
estX = zeros(N,1);
estX(estSupportSet) =  A(:,estSupportSet)\b;
end

%%%%%%%%%%%%%%%%%%%%%%%%% End of function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [estX, estSupportSet] = OMPwithInitSupport(A, b, K, initSupportSet)
%  Orthogonal Matching Pursuit Algorithm Implementation which supports partial 
% information about support-set
%

% Inputs
%   A           	: Sensing Matrix (size MxN)
%   b           	: Measurement Vector (size Mx1) b = Ax
%   K           	: sparsity level
%   initSupportSet 	: initial Support-set
% Outputs
%   estSupprtSet      	:  estimated support of x
%	estX            :  estimated sparse signal
% Purpose               : Implement Orhtogonal Matching Pursuit Algorithm 
% Reference 		: 
% Author            	: Sooraj K. Ambat
% Email              	: sooraj@ece.iisc.ernet.in/sooraj.k.ambat@gmail.com
% Address          	: Ph.D. Scholar,
%                         Statistical Signal Processing Lab,
%                         Indian Institute of Science, Bangalore, India-560 012.
% Date               	: 20 Aug 2011
% Modified 	  	: 26 May 2012

% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Short Disclaimer: this script is for educational purpose only.


%update residue with the initSupport
rk = b-A(:,initSupportSet)*(A(:,initSupportSet)\b);
estSupportSet = initSupportSet;
for i=length(initSupportSet)+1:K
    %find the atom which gives maximum correlation with current residue
    matchFilter = abs(A'*rk);
    [~, index] = max(matchFilter);
    estSupportSet = union(estSupportSet, index);
    %update residue
     rk = b-A(:,estSupportSet)*(A(:,estSupportSet)\b);
end
N = size(A,2);
estX = zeros(N,1);
estX(estSupportSet)= A(:,estSupportSet)\b;
end

%%%%%%%%%%%%%%%%%%%%%%%%% End of Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

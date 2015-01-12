% This function provides an easy way to adjust all the simulation parameters
% Variable Descriptions
%   T              	: Number of times 'x' and 'A' are independently generated
%   M              	: Number of measurements
%   L              	: Look ahead parameter
%   N              	: Ambient dimension of the signal
%   K           	: Sparsity level
%   initSupportSet 	: Initial Support-set
%   asce_omp        : Average Support Cardinality Error obtained by using OMP algorithm
%   asce_sp         : Average Support Cardinality Error obtained by using SP algorithm
%   asce_rlaomp     : Average Support Cardinality Error obtained by using Reduced LAOMP algorithm
%   srer_num		: Square of L2 norm of the signal 
%   srer_denom_omp  : Square of L2 norm of the error obtained when using OMP
%   srer_denom_sp   : Square of L2 norm of the error obtained when using SP
%   srer_denom_rlaomp: Square of L2 norm of the error obtained when using Reduced LAOMP

% Author            : Prateek B S
% Email             : bs.prateek@gmail.com
% Affiliation       : Undergrad student,
%                     Dept. of Electronics and Communication Engineering,
%                     R.V. College of Engineering, Bangalore, India.
% Info.             : This work was carried out at the Statistical Signal 
%                     Processing Lab at the Indian Institute of Science
% Date              : January 22 2014
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer: This script is for educational purposes only.

 
function [N,K,T,M,meandist,asce_omp,asce_sp,asce_laomp,asce_rlaomp,initSupportSet,L,srer_denom_omp,srer_denom_sp,srer_denom_laomp,srer_denom_rlaomp,srer_num,t_laomp,t_rlaomp]=simulation_setup()
N=500;
K=20;
T=500;
% Set the parameter M as follows to reproduce the exact same plots
% M= 30:5:95 for Gaussian Sparse Signals and Clean Measurements
% M= 30:5:140 for for Gaussian Sparse Signals and Noisy Measurements
% M= 30:5:120 for Rademacher Sparse Signals and Clean Measurements
% M= 30:5:140 for Rademacher Sparse Signals and Noisy Measurements
M=30:5:95;
meandist=0;
asce_omp=zeros(length(M),T);
asce_sp=zeros(length(M),T);
asce_laomp=zeros(length(M),T);
asce_rlaomp=zeros(length(M),T);
initSupportSet=[];
L=K/2;
srer_num=zeros(length(M),T);
srer_denom_omp=zeros(length(M),T);
srer_denom_sp=zeros(length(M),T);
srer_denom_laomp=zeros(length(M),T);
srer_denom_rlaomp=zeros(length(M),T);

t_laomp=zeros(length(M),T);
t_rlaomp=zeros(length(M),T);

end
 


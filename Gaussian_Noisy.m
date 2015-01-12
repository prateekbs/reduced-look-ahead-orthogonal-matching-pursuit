% Fig. 3 a) and Fig. 4 a) in the 'Reduced Look Ahead Orthogonal Matching
% Pursuit' paper correspond to the 2 plots generated by running this piece
% of code
% Variable Descriptions
%   x			    : Signal to be acquired (N x 1)
%   A           	: Sensing Matrix (size M x N)
%   b           	: Measurement Vector (size M x 1) b = Ax
%   T              	: Number of times 'x' and 'A' are independently generated
%   M              	: Number of measurements
%   L              	: Look ahead parameter
%   N              	: Ambient dimension of the signal
%   K           	: Sparsity level
%   asce_omp        : Average Support Cardinality Error obtained by using OMP algorithm
%   asce_sp         : Average Support Cardinality Error obtained by using SP algorithm
%   asce_rlaomp     : Average Support Cardinality Error obtained by using Reduced LAOMP algorithm
%   asce_laomp      : Average Support Cardinality Error obtained by using LAOMP algorithm
%   g               : Used to iterate for different values of M (measured vector dimension)
%   srer_omp		: Signal to Reconstruction Error Ratio (in dB) obtained by using OMP algorithm
%   srer_sp		    : Signal to Reconstruction Error Ratio (in dB) obtained by using SP algorithm
%   srer_rlaomp		: Signal to Reconstruction Error Ratio (in dB) obtained by using Reduced LAOMP algorithm
%   srer_laomp		: Signal to Reconstruction Error Ratio (in dB) obtained by using LAOMP algorithm
%   meandist        : Mean of the array of random numbers generated while creating the sensing matrix	
%   actualsupport   : Indices corresponding to the True support of the signal 
%   noisevar		: Variance of the noise added
%   initSupportSet 	: Initial Support-set
%   estSupportSetomp : estimated support of x obtained by using OMP algorithm
%   estXomp         : estimated sparse signal obtained by using OMP algorithm
%   estSupportSetsp  : estimated support of x obtained by using SP algorithm
%   estXsp          : estimated sparse signal obtained by using SP algorithm
%   estSupportSetrlaomp: estimated support of x obtained by using Reduced LAOMP algorithm
%   estSupportSetlaomp: estimated support of x obtained by using Reduced LAOMP algorithm
%   estXrlaomp      : estimated sparse signal obtained by using Reduced LAOMP algorithm
%   estXlaomp      : estimated sparse signal obtained by using Reduced LAOMP algorithm
%   Union_set		: Union of the support set estimated by OMP and SP 
%   Ik			    : Intersection of the support set estimated by OMP and SP
%   srer_num		: Square of L2 norm of the signal 
%   srer_denom_omp  : Square of L2 norm of the error obtained when using OMP
%   srer_denom_sp   : Square of L2 norm of the error obtained when using SP
%   srer_denom_rlaomp: Square of L2 norm of the error obtained when using Reduced LAOMP
%   srer_denom_laomp: Square of L2 norm of the error obtained when using LAOMP

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

clear all; close all; clc;
% The function simulation_setup() sets up all the simulation parameters. To
% see what actual values are set, please open the simulation_setup.m file
[N,K,T,M,meandist,asce_omp,asce_sp,asce_laomp,asce_rlaomp,initSupportSet,L,srer_denom_omp,srer_denom_sp,srer_denom_laomp,srer_denom_rlaomp,srer_num,t_laomp,t_rlaomp]=simulation_setup();
% Used to set the seed to the default value. A very helpful tool for
% reproducible research
RandStream.setDefaultStream(RandStream('mt19937ar','seed',0));
for g=1:length(M)
    for t=1:T
        % Generate Measurement matrix and normalize all the columns
        A=normc(normrnd(meandist,1/M(g),M(g),N));
        % Generate Signal and its corresponsing support. In the first argument of the function, Insert
        % 'BernoulliRndSign' for Rademacher Sparse signals and 
        % 'StdGaussian' for Gaussian Sparse Signals.
        [ x, actualsupport] = generateSignal('StdGaussian', N, K  );
        y=A*x;
        % Insert noise to ensure signal to measurement noise (SMNR) of 15 dB
        noiseVar = (K/M(g))*10^(-15/10);
        y=y+sqrt(noiseVar)*randn(M(g),1);
        
        temp=cputime;
        [estXlaomp,estSupportSetlaomp]=laomp(A,y,K,L);
        % Store CPU time for LAOMP
        t_laomp(g,t)=cputime-temp;
        temp=cputime; 
        % Use functions to compute the estimated value of signal using OMP
        %, SP and LAOMP.
        [estXomp, estSupportSetomp] = OMPwithInitSupport(A, y, K, initSupportSet);
        [estXsp, estSupportSetsp] = SPwithInitialSupport(A, y, K, initSupportSet);
        % Find the union and intersection of the support sets
        Union_set=union(estSupportSetomp,estSupportSetsp);
        Ik=intersect(estSupportSetomp,estSupportSetsp);
        % Use a Reduced LAOMP algorithm to search in only in the Union
        % space estimated by OMP and SP. Obtain the estimated value of
        % signal using the proposed algorithm
        [estXrlaomp,estSupportSetrlaomp]=modified_laomp(A,y,K,L,Ik,Union_set);
        % Store CPU time for RLAOMP
        t_rlaomp(g,t)=cputime-temp;
        % Compute the ASCE and SRER for OMP,SP,LAOMP and the proposed Reduced LAOMP for
        % this iteration and store in a matrix.
        asce_omp(g,t)=1-(1/K)*(length(intersect(actualsupport,estSupportSetomp)));
        asce_sp(g,t)=1-(1/K)*(length(intersect(actualsupport,estSupportSetsp)));
        asce_laomp(g,t)=1-(1/K)*(length(intersect(actualsupport,estSupportSetlaomp)));
        asce_rlaomp(g,t)=1-(1/K)*(length(intersect(actualsupport,estSupportSetrlaomp)));
        srer_num(g,t)=power(norm(x),2);
        srer_denom_omp(g,t)=power(norm(x-estXomp),2);
        srer_denom_sp(g,t)=power(norm(x-estXsp),2);
        srer_denom_laomp(g,t)=power(norm(x-estXlaomp),2);
        srer_denom_rlaomp(g,t)=power(norm(x-estXrlaomp),2);
    end
    
end
% Average over all iterations at a particular fraction of measurement value
asce_omp=mean(asce_omp.');
asce_sp=mean(asce_sp.');
asce_laomp=mean(asce_laomp.');
asce_rlaomp=mean(asce_rlaomp.');


% Plot fraction of measurements v/s Average Support Cardinality Error
plot(M./N,asce_omp,'-*b',M./N,asce_sp,'-+r',M./N,asce_laomp,'-.magenta',M./N,asce_rlaomp,'-oblack');
legend('OMP ','SP ','LAOMP','Reduced LAOMP ');
xlabel(' Fraction of Measurements');
ylabel('Average Support Cardinality Error');
title('ASCE for Gaussian Sparse Signals (SMNR= 15dB)');
grid;
axis tight;
figure;
% Average over all iterations at a particular fraction of measurement value
srer_num=mean(srer_num.');
srer_denom_omp1=mean(srer_denom_omp.');
srer_denom_sp1=mean(srer_denom_sp.');
srer_denom_laomp1=mean(srer_denom_laomp.');
srer_denom_rlaomp1=mean(srer_denom_rlaomp.');
srer_omp=10*log10(srer_num)-10*log10(srer_denom_omp1);
srer_sp=10*log10(srer_num)-10*log10(srer_denom_sp1);
srer_laomp=10*log10(srer_num)-10*log10(srer_denom_laomp1);
srer_rlaomp=10*log10(srer_num)-10*log10(srer_denom_rlaomp1);

% Plot fraction of measurements v/s Signal to Reconstruction Error Ratio
plot(M./N,srer_omp,'-*b',M./N,srer_sp,'-+r',M./N,srer_laomp,'-.magenta',M./N,srer_rlaomp,'-oblack');
legend('OMP ','SP ','LAOMP','Reduced LAOMP ');
xlabel(' Fraction of measurements');
ylabel('Signal to Reconstruction Error Raio (SRER) in dB');
title('SRER for Gaussian Sparse signals (SMNR= 15dB)');
grid;
axis tight;

% Average over all iterations at a particular fraction of measurement value
t_laomp=mean(t_laomp.'); 
t_rlaomp=mean(t_rlaomp.'); 

save('Gaussian Noisy');

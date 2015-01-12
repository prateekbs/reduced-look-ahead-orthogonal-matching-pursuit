% This script executes the residue fucntion described in the 'Look Ahead Orthogonal Matching Pursuit'. For each individual choice of L, the   % algorithm is executed like a normal OMP until the sparsity level is met. The best atom among L atoms that results in the minimum norm of 
% fitting residual is then chosen. 

% Variable Descriptions
%   A           	: Sensing Matrix (size M x N)
%   y           	: Measurement Vector (size M x 1) y = Ax
%   K           	: Sparsity level
%   rk                  : Final residue
%   I 			: Initial estimated support set
%   i			: One of the L atoms picked in the LAOMP algorithm

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
function rk=laomp_residue(A,y,K,I,i)

    Ik=union(I,i);
    estX=zeros(size(A,2),1);
    estX(Ik)=A(:,Ik)\y;
    rk=y-A*estX;
    k=length(Ik);
   for k=k:K
	% Finding index of the atom from 'A' that is maximally correlated to residue
     [~,ik]=max(abs(A'*rk));
     Ik=union(ik,Ik);
     estX(Ik)=A(:,Ik)\y;
     rk=y-A*estX;  
   end
    
end

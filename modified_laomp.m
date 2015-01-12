% This function is a modified version of LAOMP algorithm proposed in the paper 'Look Ahead Orthogonal Matching Pursuit'. 
% The modified LAOMP does not search the entire space spanned by the columns of A for K atoms. It searches for K-Ik elements in the space
% spanned by the Union_Set only.
% Variable Descriptions
%   A           	: Sensing Matrix (size M x N)
%   y           	: Measurement Vector (size M x 1) y = Ax
%   L              	: Look ahead parameter
%   K           	: Sparsity level
%   estX                : estimated support of x obtained by using Reduced LAOMP algorithm
%   estSupportSetlaomp  : estimated support of x obtained by using Reduced LAOMP algorithm
%   Union_set		: Union of the support set estimated by OMP and SP 
%   Ik			: Intersection of the support set estimated by OMP and SP

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

function [estX,estSupportSetlaomp]=modified_laomp(A,y,K,L,Ik,Union_set)
    n=zeros(L,1);
    estX=zeros(size(A,2),1);
    estX(Ik)=A(:,Ik)\y;
    rk=y-A*estX;
    k=length(Ik);
    % Search only for 'k - length of intersection' atoms	   
	for k=k:K
        [~,indices]=sort(abs(A(:,Union_set)'*rk),'descend');
        indices=Union_set(indices);
                
	j=indices(1:L);
	% Run the residue calculation function for each of the L atoms        
	for l=1:L
            rr=laomp_residue(A,y,K,Ik,j(l));
            n(l)=norm(rr,2);
        end
        [~,indices]=sort(n,'ascend');
	% Choose the atom that results in the lowest residue
        l=indices(1);
        ik=j(l);
	% Add chosen atom to the estimated support set
        Ik=union(Ik,ik);
        estX(Ik)=A(:,Ik)\y;
        rk=y-A*estX;
        
    end
   estSupportSetlaomp=find(estX);
end

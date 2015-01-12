% This fucntion is a modified version of LAOMP proposed by Saikat
% Chatterjee and group. The modified LAOMP does not search the entire space
% spanned by the columns of A. It searches for K-Ik elements in the space
% spanned by the Union_Set only.
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

function [estX,estSupportSetlaomp]=laomp(A,y,K,L)
    n=zeros(L,1);
    estX=A\y;
    rk=y-A*estX;
    Ik=[];
    for k=1:K
        [~,indices]=sort(abs(A'*rk),'descend');
        j=indices(1:L);
        for l=1:L
            rr=laomp_residue(A,y,K,Ik,j(l));
            n(l)=norm(rr,2);
        end
        [~,indices]=sort(n,'ascend');
        l=indices(1);
        ik=j(l);
        Ik=union(Ik,ik);
        estX=zeros(size(A,2),1);
        estX(Ik)=A(:,Ik)\y;
        rk=y-A*estX;
        
    end
   estSupportSetlaomp=find(estX);
end

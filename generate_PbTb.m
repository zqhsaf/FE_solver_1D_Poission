function [Pb,Tb]=generate_PbTb(left,right,number_of_element,basis_type)

N=number_of_element;

if basis_type==101

    Pb=left:1/N:right;

    Tb=[1:N;2:N+1];

elseif basis_type==102

    Pb=left:1/(2*N):right;

    Tb=[1:2:2*N-1;3:2:2*N+1;2:2:2*N];
        
end
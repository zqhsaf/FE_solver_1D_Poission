function [Pb,Tb]=generate_PbTb(left,right,number_of_element,basis_type)

N=number_of_element;

if basis_type==101

    Pb=left:1/N:right;

    Tb=[1:N;2:N+1];
        
end
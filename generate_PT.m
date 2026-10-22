function [P,T]=generate_PT(left,right,N)

    P=left:1/N:right;

    T=[1:N;2:N+1];

end
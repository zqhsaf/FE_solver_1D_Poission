function [A,b] = treat_Dirichlet_boundary(A, b, boundarynodes, Pb, g, dirichlet_tag)

if nargin < 6 || isempty(dirichlet_tag)

    dirichlet_tag = 1;     

end

nbn = size(boundarynodes, 2);

for k = 1 : nbn

    if boundarynodes(1, k) == dirichlet_tag

        i = boundarynodes(2, k);     

        A(i, :) = 0;                   

        A(i, i) = 1;                   

        b(i) = g(i);

    end

end

end

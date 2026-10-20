function [Gauss_weights,Gauss_nodes]=generate_Gauss_local(vertices,Gpn)

if Gpn < 1

    error('n>=1'); 

end

h=vertices(2)-vertices(1);

k = (1:Gpn-1)';

b = k ./ sqrt(4*k.^2 - 1);

T = diag(b,1) + diag(b,-1);

[V,D] = eig(T);

[~,idx] = sort(diag(D));

Gauss_nodes_origin = diag(D(idx,idx));

Gauss_weights_origin = 2 * V(1,idx).^2;

Gauss_nodes  = h/2 * Gauss_nodes_origin + (vertices(1)+vertices(2))/2;

Gauss_weights= h/2 * Gauss_weights_origin;

end
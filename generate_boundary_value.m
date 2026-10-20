function g=generate_boundary_value(Pb,boundarynodes,value1,value2)

g=sparse(1,length(Pb));

g(boundarynodes(2,1))=value1;

g(boundarynodes(2,2))=value2;

end
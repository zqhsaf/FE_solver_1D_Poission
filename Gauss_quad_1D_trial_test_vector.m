function int_value=Gauss_quad_1D_trial_test_vector(f_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_test,basis_index_test,basis_der_x_test)

Gpn=5; 
int_value=0;

for k=1:Gpn
    
    int_value=int_value+Gauss_weights(k)*feval(f_fun,Gauss_nodes(k))*FE_basis_local_fun_1D(Gauss_nodes(k),vertices,basis_type_test,basis_index_test,basis_der_x_test);

end
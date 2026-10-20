function b=assemble_vector_1D(f_fun,vector_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,basis_der_x_test)

b=zeros(vector_size);

Gpn=5;

for n=1:number_of_elements

    vertices=P(:, T(:,n));   
    [Gauss_weights,Gauss_nodes]=generate_Gauss_local(vertices,Gpn);

        for beta=1:number_of_local_basis_fun_test

            int_value=Gauss_quad_1D_trial_test_vector(f_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_test,beta,basis_der_x_test);

            b(Tb_test(beta,n),1)=b(Tb_test(beta,n),1)+int_value;

        end

end
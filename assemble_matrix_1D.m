function A=assemble_matrix_1D(coe_fun,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test)

A=sparse(matrix_size(1),matrix_size(2));

for n=1:number_of_elements

    Gpn=5;
    vertices=P(:, T(:,n));   
    [Gauss_weights,Gauss_nodes]=generate_Gauss_local(vertices,Gpn);

    for alpha=1:number_of_local_basis_fun_trial

        for beta=1:number_of_local_basis_fun_test

            int_value=Gauss_quad_1D_trial_test_matrix(coe_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_trial,alpha,basis_der_x_trial,basis_type_test,beta,basis_der_x_test);

            %S(beta,alpha)=int_value; 
            i=Tb_test(beta,n);
            j=Tb_trial(alpha,n);
            A(i,j)=A(i,j)+int_value;
        
        end
    end
end
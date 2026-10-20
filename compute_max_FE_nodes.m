function max_error = compute_max_FE_nodes(exact_fun, Pb, uh)


u_exact = feval(exact_fun, Pb);

max_error = max(abs(uh(:) - u_exact(:)));

end

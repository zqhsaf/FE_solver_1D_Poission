%% Robin(第三类)边界例题: u(0)=0 (Dirichlet), u'(1)+q_b*u(1)=p_b (Robin)
%  -d/dx(e^x du/dx) = -e^x[cos(x)-2sin(x)-xcos(x)-xsin(x)], 0<=x<=1
%  精确解 u(x)=x cos(x)
%
%  Robin 边界处理(弱形式):
%    ∫ c u'v' dx + q_b c(b) u(b)v(b) = ∫ f v dx + p_b c(b) v(b)
%  即: 刚度矩阵 A(b,b) 加 q_b*c(b)，载荷向量 b(b) 加 p_b*c(b)
%
%  为用精确解 u=x cos(x) 验证，取 q_b=1，则
%    p_b = u'(1) + q_b*u(1) = (cos1-sin1) + cos1 = 2cos1-sin1
clear; clc;

qb=1;
pb=2*cos(1)-sin(1);
N=4;

%% 1. 一次元(101)求解
u1=FE_solver_1D_Poission(0,1,101,101,N,0,pb,'function_u_exact',1,3,0,qb);

Pb1=0:1/N:1;
fprintf('\n===== 一次元 N=%d 节点结果 =====\n',N);
disp(table(Pb1',u1,Pb1'.*cos(Pb1'),u1-Pb1'.*cos(Pb1'), ...
    'VariableNames',{'x','u_h','u_exact','error'}));

%% 2. 二次元(102)求解
u2=FE_solver_1D_Poission(0,1,102,102,N,0,pb,'function_u_exact',1,3,0,qb);

Pb2=0:1/(2*N):1;
fprintf('===== 二次元 N=%d 节点结果 =====\n',N);
disp(table(Pb2',u2,Pb2'.*cos(Pb2'),u2-Pb2'.*cos(Pb2'), ...
    'VariableNames',{'x','u_h','u_exact','error'}));

%% 3. 直接检验 Robin 条件: u_h'(1)+q_b*u_h(1) 应逼近 p_b
[P,T]=generate_PT(0,1,N);
[~,Tb1]=generate_PbTb(0,1,N,101);
[~,Tb2]=generate_PbTb(0,1,N,102);
vertices=P(:,T(:,N));
rb1=local_FE_function_1D(1,u1(Tb1(:,N)),vertices,101,1)+qb*u1(end);
rb2=local_FE_function_1D(1,u2(Tb2(:,N)),vertices,102,1)+qb*u2(end);
fprintf('Robin 条件检验: u_h''(1)+q_b*u_h(1) 一次元=%.8f, 二次元=%.8f, p_b=%.8f\n', ...
    rb1,rb2,pb);

%% 4. 一致性检验: q_b=0 的 Robin 应与 Neumann 结果完全相同
uR=FE_solver_1D_Poission(0,1,102,102,N,0,cos(1)-sin(1),[],1,3,0,0);
uN=FE_solver_1D_Poission(0,1,102,102,N,0,cos(1)-sin(1),[],1,2);
fprintf('q=0 的 Robin 与 Neumann 结果之差 = %.3e (应为 0)\n',norm(uR-uN));

%% 5. 收敛阶验证
fprintf('\n===== 一次元收敛阶 (理论: 节点/L2 为 O(h^2), H1 为 O(h)) =====\n');
print_order(101,qb,pb);
fprintf('\n===== 二次元收敛阶 (理论: L2 为 O(h^3), H1 为 O(h^2)) =====\n');
print_order(102,qb,pb);

function print_order(basis_type,qb,pb)

fprintf('%4s %12s %8s %12s %8s %12s %8s\n', ...
    'N','max节点误差','阶','L2误差','阶','H1误差','阶');

em_prev=0; eL_prev=0; eH_prev=0;

for N=[4 8 16 32]

    [P,T]=generate_PT(0,1,N);
    [Pb,Tb]=generate_PbTb(0,1,N,basis_type);

    u=FE_solver_1D_Poission(0,1,basis_type,basis_type,N,0,pb,[],1,3,0,qb);

    em=compute_max_FE_nodes('function_u_exact',Pb,u);
    eL=compute_Hs_error('analytic_solution',0,N,P,T,Tb,u,basis_type);
    eH=compute_Hs_error('analytic_solution',1,N,P,T,Tb,u,basis_type);

    if em_prev==0
        fprintf('%4d %12.3e %8s %12.3e %8s %12.3e %8s\n',N,em,'-',eL,'-',eH,'-');
    else
        fprintf('%4d %12.3e %8.2f %12.3e %8.2f %12.3e %8.2f\n', ...
            N,em,log2(em_prev/em),eL,log2(eL_prev/eL),eH,log2(eH_prev/eH));
    end

    em_prev=em; eL_prev=eL; eH_prev=eH;
end

end

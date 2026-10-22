function [A,b]=treat_Robin_boundary_1D(coe_fun,A,b,boundarynodes,Pb,boundaryvalues,boundaryq)
% treat_Robin_boundary_1D 处理 Robin(第三类)边界 u'(x)+q(x)u(x)=p(x)
%
%   弱形式经分部积分并代入边界条件:
%     右端点 b: u'(b)=p_b-q_b u(b)
%       q_b c(b) u(b)v(b) + ∫ c u'v' dx = ∫ f v dx + p_b c(b) v(b)
%       => A(b,b) += +q_b c(b),  b(b) += +p_b c(b)
%     左端点 a: u'(a)=p_a-q_a u(a)
%       ∫ c u'v' dx - q_a c(a) u(a)v(a) = ∫ f v dx - p_a c(a) v(a)
%       => A(a,a) += -q_a c(a),  b(a) += -p_a c(a)
%
%   输入:
%       coe_fun        : 系数函数 c(x)
%       A,b            : 已装配的刚度矩阵、载荷向量
%       boundarynodes  : 2×nbn 边界节点信息矩阵，第1行为边界类型(1=Dirichlet,2=Neumann,3=Robin)
%       Pb             : 节点坐标
%       boundaryvalues : 按节点编号存放的 p 值（边界值）
%       boundaryq      : 按节点编号存放的 q 值（Robin 系数）
%
%   输出:
%       A,b : 加入 Robin 边界项后的刚度矩阵、载荷向量

nbn=size(boundarynodes,2);

for k=1:nbn

    if boundarynodes(1,k)==3

        i=boundarynodes(2,k);

        if k==1
            sign_coef=-1;   % 左端点 a: 刚度/载荷项均带负号
        else
            sign_coef=1;    % 右端点 b: 刚度/载荷项均带正号
        end

        A(i,i)=A(i,i)+sign_coef*feval(coe_fun,Pb(i))*boundaryq(i);

        b(i,1)=b(i,1)+sign_coef*feval(coe_fun,Pb(i))*boundaryvalues(i);

    end

end

end

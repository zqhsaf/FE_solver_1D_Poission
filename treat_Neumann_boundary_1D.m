function b=treat_Neumann_boundary_1D(coe_fun,b,boundarynodes,Pb,boundaryvalues)
% treat_Neumann_boundary_1D 处理 Neumann 边界 u'(x)=r
%
%   弱形式经分部积分后:
%       c(a)u'(a)v(a) - c(b)u'(b)v(b) + ∫ c u'v' dx = ∫ f v dx
%   移项得:
%       ∫ c u'v' dx = ∫ f v dx - r_a c(a)v(a) + r_b c(b)v(b)
%   故:
%       左端点(a)向载荷向量减去 r_a c(a)，右端点(b)向载荷向量加上 r_b c(b)
%
%   输入:
%       coe_fun        : 系数函数 c(x)
%       b              : 已装配的载荷向量
%       boundarynodes  : 2×nbn 边界节点信息矩阵，第1行为边界类型(1=Dirichlet,2=Neumann)
%       Pb             : 节点坐标
%       boundaryvalues : 按节点编号存放的边界值向量
%
%   输出:
%       b : 加入 Neumann 边界项后的载荷向量

nbn=size(boundarynodes,2);

for k=1:nbn

    if boundarynodes(1,k)==2

        i=boundarynodes(2,k);

        if k==1
            sign_coef=-1;   % 左端点 a: -r_a c(a)v(a)
        else
            sign_coef=1;    % 右端点 b: +r_b c(b)v(b)
        end

        b(i,1)=b(i,1)+sign_coef*feval(coe_fun,Pb(i))*boundaryvalues(i);

    end

end

end

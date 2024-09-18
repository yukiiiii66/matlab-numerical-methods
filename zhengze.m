function [x,A,r] = zhengze(m,n,b)
%20211081037 2103 俞昊然
%此为正则化方法求解最小二乘问题的程序
%m为系数矩阵的行数，n为列数，b为常向量，b的维数应等于m
if m<n
    error('你所输入的问题一定有解，不适合寻找最小二乘解');
end
%随机生成矩阵A
A=randn(m,n);
r=rank(A);
if r<n
    fprintf('你所输入的矩阵为秩亏损型')'
else
    A1=A'*A;
    B=A'*b;
    %使用cholesky分解进行计算
    [x] = cholesky_auto(A1,B);
    %计算残量
    r=norm((A'*(A*x-b)));
end


end
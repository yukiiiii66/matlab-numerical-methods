function [namb,x] = fmifa(A,u,eps,k)
%20211081037 俞昊然 数计2103
%返回namb为A模最小的特征值，x为其对应的特征向量
%A为求特征值的矩阵，u为初始向量，eps为收敛精度，k为最大迭代步数
e=1;
iter=0;
%进行反幂法迭代
while e>eps&&iter<k
    [~,index]=max(abs(u));
    y=u./u(index);
    %vk=A逆*zk-1求vk
    u1=A\y;
    e=norm(u1-u);
    %更新迭代
    u=u1;
    iter=iter+1;
end
[~,index]=max(abs(u));
%输出模最小的特征值
namb=1/u(index);
%输出特征向量
x=y;
end


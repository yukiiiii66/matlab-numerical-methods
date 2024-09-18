function [x] = SOR2(w,eps,A)
%该函数为带松弛因子的SOR迭代方法，w为松弛因子，eps为终止准则的检验精度
%20211081037-俞昊然
%A为第二题用数据做出的矩阵
s1=size(A);
s1=s1(1);
%设定初始点
b=A*zeros(s1,1);
x0=ones(s1,1);
x=x0;
iter=0;normg=norm(A*x-b)
while normg>eps && iter<40000
iter=iter+1; 

    for i=1:s1
        if A(i,i)==0
            x(i)=0;
        else
           s=0;d=0;
           %考虑在当前迭代次数中，已生成的新迭代点
           if i>1
                s=A(i,1:i-1)*x(1:i-1);
           else
           end
        %考虑在当前迭代次数中，还未生成的迭代点
           if i<s1
                d=A(i,i+1:end)*x(i+1:end);
           else
           end
          %计算当前迭代中，第i个点的新迭代值
           x(i)=(1-w)*x(i)+w/A(i,i)*(b(i)-s-d);
        end
    end
    %计算当前迭代点残量
    normg=norm(A*x-b);
    fprintf('当前迭代次数为%d 当前迭代点的参量为 %d\n',iter,normg)
end
end


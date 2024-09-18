function [v,beta] = Householder1(x)
%20211081037 2103 俞昊然
%对一个列向量进行househodler变换 使得后n-1个元素全部变为0 Hx=a e1
gama=norm(x,inf);
v=x/gama;
alpha=norm(v);
new=v(1)+sign(v(1))*alpha;
beta=new^2/((alpha+abs(v(1)))*alpha);
v(1)=1;v(2:length(x))=v(2:length(x))/new;
end


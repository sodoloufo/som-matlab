function y = sumignorenan(x,z)

%%%%%%%%%%%%%%%%%% Daouda Diouf %%%%%%%%%%%%%%%%%%%%%%
%Somme de vecteurs en ignorant les NAN 
if size(x)~=size(z)
    error('dimension of matrix must be the same')
end

mat1=x(:);
mat2=z(:);
y=nan*ones(length(mat1),1);
for i=1:length(mat1)
    if (isnan(mat1(i)) & ~isnan(mat2(i)))
        y(i)=mat2(i);
    elseif (isnan(mat2(i)) & ~isnan(mat1(i)))
        y(i)=mat1(i);
    elseif (~isnan(mat2(i)) & ~isnan(mat1(i)))
        y(i)=mat1(i)+mat2(i);
    end
end
y=reshape(y,size(x,1),size(x,2));

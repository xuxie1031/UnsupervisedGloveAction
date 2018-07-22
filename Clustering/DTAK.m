%%%%%%%%%%%%%%%%%%%
% DTAK similarity %
%%%%%%%%%%%%%%%%%%%

function [gamma,U]=DTAK(kernel_b, gamma)
    gamma=0;
    [dimy,dimx]=size(kernel_b);
    U=zeros(dimy,dimx);
    U(1,1)=2*kernel_b(1,1);
    
    for i=1:dimy
        for j=1:dimx
            if i==1&&j==1
                continue;
            elseif i==1
                U(i,j)=U(i,j-1)+kernel_b(i,j);
            elseif j==1
                U(i,j)=U(i-1,j)+kernel_b(i,j);
            else
                U(i,j)=max(max(U(i,j-1)+kernel_b(i,j),U(i-1,j)+kernel_b(i,j)),U(i-1,j-1)+2*kernel_b(i,j));
            end
        end
    end
    gamma=single(U(dimy,dimx))/(dimy+dimx);
end


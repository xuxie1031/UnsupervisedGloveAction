%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate kernel among two features %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel_b = kernel_binary(X, Y, sigma_, threshold)
    [len_Y,~]=size(Y);
    [len_X,~]=size(X);
    
    kernel_b=zeros(len_Y,len_X);
    for i=1:len_Y
        for j=1:len_X
            kernel_b(i,j)=exp(-norm(Y(i,:)-X(j,:))/(2*sigma_^2))>threshold;
        end
    end
end


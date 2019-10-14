function x_present=input_generator(x,order,n)
    x_present=zeros(order+1,1);
    for j=1:order+1,
        if((n-j) >= 0) 
            x_present(j,1) = x(n-j+1);
        else
            x_present(j,1) =0;                      
        end
    end
end
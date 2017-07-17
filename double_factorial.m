function x = double_factorial(n)

    if mod(n,2) == 0
       
        x = 1;
        for k = 1:(n/2)
           x = x*(2*k); 
        end
        
    else
        
        x = 1;
        for k = 1:(n+1)/2
           x = x*(2*k-1); 
        end        
        
    end

end
function mu = piecewise_linear_shrinker(x,a,b)
%implements piecewise linear shrinkage: 
%returns a value that equals zero on (0,a), identity o (b,Inf), and linear
%interpolation in between

J = length(x);
mu = zeros(J,1);
for i=1:J
   y = x(i);
   
   if (y <=a) 
       mu(i) =0;
   else
       if (y>=b)
           mu(i) = y;
       else
           mu(i) = b/(b-a)*(y-a);
       end
   end
end
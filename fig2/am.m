function annual_mean=am(x)
% calculates annual means from monthly means
% assumptions: time is multiple of 12 and first dimension in case of
% multidimensional input matrix (e.g., a time, lat, lon field)


if ndims(x) == 3
    
    for j=1:(length(x)/12)
        annual_mean(j,:,:) = (31*x((12*j)-12+1,:,:)+28*x((12*j)-12+2,:,:)+31*x((12*j)-12+3,:,:)+30*x((12*j)-12+4,:,:)+31*x((12*j)-12+5,:,:)+30*x((12*j)-12+6,:,:)+31*x((12*j)-12+7,:,:)+31*x((12*j)-12+8,:,:)+30*x((12*j)-12+9,:,:)+31*x((12*j)-12+10,:,:)+30*x((12*j)-12+11,:,:)+31*x((12*j)-12+12,:,:))/365;
    end
    
elseif ndims(x) == 2
    
    tmp = size(x);

    if tmp(1) == 1 || tmp(2) == 1
        for j=1:(length(x)/12)
            annual_mean(j) = (31*x((12*j)-12+1)+28*x((12*j)-12+2)+31*x((12*j)-12+3)+30*x((12*j)-12+4)+31*x((12*j)-12+5)+30*x((12*j)-12+6)+31*x((12*j)-12+7)+31*x((12*j)-12+8)+30*x((12*j)-12+9)+31*x((12*j)-12+10)+30*x((12*j)-12+11)+31*x((12*j)-12+12))/365;
        end
    else    
        for j=1:(length(x)/12)
            annual_mean(j,:) = (31*x((12*j)-12+1,:)+28*x((12*j)-12+2,:)+31*x((12*j)-12+3,:)+30*x((12*j)-12+4,:)+31*x((12*j)-12+5,:)+30*x((12*j)-12+6,:)+31*x((12*j)-12+7,:)+31*x((12*j)-12+8,:)+30*x((12*j)-12+9,:)+31*x((12*j)-12+10,:)+30*x((12*j)-12+11,:)+31*x((12*j)-12+12,:))/365;
        end
    end
else
    
    'too many dimension for this function!'
    
end
    
return
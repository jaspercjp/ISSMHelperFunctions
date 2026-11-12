function [arr] = fill_nan(arr)
    for i=2:length(arr)
        if isnan(arr(i))
            arr(i)=arr(i-1);
        end
    end

    for i=length(arr)-1:-1:1
        if isnan(arr(i))
            arr(i)=arr(i+1);
        end
    end
end
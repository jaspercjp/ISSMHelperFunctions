function s = dec2str(num)
    if num<1 && num>0.1
        s = "0p" + floor(100*num);
    elseif num<1 && num<0.1
        s="0p0"+floor(100*num);
    else
        tmp1 = string(floor((floor(100*num) - mod(100*num, 100))/100 + 1e-6));
        s = tmp1 + "p" + string(floor(mod(100*num,100)));
    end
end
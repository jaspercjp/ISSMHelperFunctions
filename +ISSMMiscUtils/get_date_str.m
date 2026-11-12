function date_str = get_date_str()
    [hr,minute,~]=hms(datetime);
    if minute<10
        minute_string = "0"+string(minute);
    else
        minute_string = string(minute);
    end

    if hr<10
        hr_string = "0" + string(hr);
    else
        hr_string = string(hr);
    end
    date_str = string(yyyymmdd(datetime("today")))+"_"+hr_string+minute_string;
end

function [ST,ST2]=write_soil_rout_mvl(ST)
ST=ST/1000; %mm --> m
if ST <= 0.3
    ST = 0.3; 
end
ST2=ST;


end


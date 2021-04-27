function [ climateOS_INPUT,Pout,climateOS_Pori] = write_rain_rout365_mvl(climateOS,time,extraP4cell)
climateOS_INPUT=climateOS; %if not indexed w/ j, climateOS will result in different variables depending on the parfor iteration

climateOS_INPUT(:,6)=climateOS(:,6)'+extraP4cell;
climateOS_Pori=climateOS(:,6);

Pout=sum(climateOS_INPUT(:,6));


end




function [ST] = get_soil_rout(time,ii,jj,xx1,yy1,sedexport,marshmask,testmodeSTshort,datapath)

path = [datapath 'ST/'];
% file = sprintf('%s%d%s','Soil Thickness in mm',(fix(time/10)*10)+10,'_1_1'); % To get same ST map for each 10 yearly period!! you need to do 4000runs to get all the years done, or take 10 yearly climate averages. . . 
file = sprintf('%s%d%s','Soil Thickness in mm',time*10,'_1_1');
if testmodeSTshort == 1
    file = sprintf('%s%d%s','Soil Thickness in mm',time*200,'_1_1');
end
if testmodeSTshort == 2
    if time ==1
        stmaparr=[repmat(10,1,10),repmat(20,1,10)];
        file = sprintf('%s%d%s','Soil Thickness in mm',stmaparr(time),'_1_1');
    else
        stmaparr=[repmat(10,1,10),repmat(20,1,10)];
        file = sprintf('%s%d%s','Soil Thickness in mm',stmaparr(time-1),'_1_1');
    end
end

fname= [path file];
[x,row,col,hres,vres] = Read_Idrisi(fname);
% flag = 0;
% if mask(i+xx1,j+yy1) == 0   %indien full GGZ run, moet dit ook aangepast worden naar +xx1 +yy1 %% niet meer nodig gezien landuse al gecheckt wordt in vorige functie.
%     flag = 1;
% end
% if flag == 0; % we want to know the CY, so read ST
    ST=x(jj+yy1-1,ii+xx1-1); %xx1 and yy1 transform i and j into real actual coordinates for the gravgaz ST map. 
%     if ST == 0
%         ST = 0.1;% Aq doesn't work when using a ST of zero.
%     end


%  if marshmask(jj+yy1-1,ii+xx1-1)==1   %% only 4 GGZ
%     ST=1500+(((sedexport(time)/1.35)/(numel(marshmask(marshmask==1))*400))*1000); %sedexport =t/j, -> m³/m² -> m->mm
%  end
%  
 
% else
%     ST=0;%we don't want to know the CY, so we'll put the ST on zero so that CY is also zero. Not flexible, TO DO !
% end

%niet nodig...
% ST=round(ST/100)*100; % voor de soil comparments niet dikker te maken dan de soil zelf is, testje...
end


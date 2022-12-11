%%-------------------------------------------------------------------------
% Usage: the Chauvenet's Criterion method to identify outliers
%%-------------------------------------------------------------------------
% Substrate: MUF-P
%%-------------------------------------------------------------------------
clear 
close all 
%%
% read data 
lat = xlsread('APA_data.xlsx',1,'B:B');
lon = xlsread('APA_data.xlsx',1,'C:C');
APA_bulk = xlsread('APA_data.xlsx',1,'J:J');
APA_dissolved = xlsread('APA_data.xlsx',1,'K:K');
APA_particulate = xlsread('APA_data.xlsx',1,'L:L');
APA_bacterial = xlsread('APA_data.xlsx',1,'M:M');
APA_phytoplankton = xlsread('APA_data.xlsx',1,'N:N');
APA_trichodesmium = xlsread('APA_data.xlsx',1,'O:O');
% Set the data -999 to the null value NaN
APA_bulk(APA_bulk==-999)=NaN;
APA_dissolved(APA_dissolved==-999)=NaN;
APA_particulate(APA_particulate==-999)=NaN;
APA_bacterial(APA_bacterial==-999)=NaN;
APA_phytoplankton(APA_phytoplankton==-999)=NaN;
APA_trichodesmium(APA_trichodesmium==-999)=NaN;
%%
APA_bulk(APA_bulk==0)= NaN;
APA_dissolved(APA_dissolved==0)= NaN;
APA_particulate(APA_particulate==0)= NaN;
APA_bacterial(APA_bacterial==0)= NaN;
APA_phytoplankton(APA_phytoplankton==0)= NaN;
APA_trichodesmium(APA_trichodesmium==0)= NaN;
% log transform
APA_bulk_log=log10(APA_bulk);
APA_dissolved_log=log10(APA_dissolved);
APA_particulate_log=log10(APA_particulate);
APA_bacterial_log=log10(APA_bacterial);
APA_phytoplankton_log=log10(APA_phytoplankton);
APA_trichodesmium_log=log10(APA_trichodesmium);
% calculate mean value
APA_bulk_log_m=mean(APA_bulk_log(:,1),'omitnan');
APA_dissolved_log_m=mean(APA_dissolved_log(:,1),'omitnan');
APA_particulate_log_m=mean(APA_particulate_log(:,1),'omitnan');
APA_bacterial_log_m=mean(APA_bacterial_log(:,1),'omitnan');
APA_phytoplankton_log_m=mean(APA_phytoplankton_log(:,1),'omitnan');
APA_trichodesmium_log_m=mean(APA_trichodesmium_log(:,1),'omitnan');
APA_log_mean=[APA_bulk_log_m,APA_dissolved_log_m,APA_particulate_log_m,APA_bacterial_log_m,APA_phytoplankton_log_m,APA_trichodesmium_log_m];
% calculate standard deviation
APA_bulk_log_std=std(APA_bulk_log(:,1),'omitnan');
APA_dissolved_log_std=std(APA_dissolved_log(:,1),'omitnan');
APA_particulate_log_std=std(APA_particulate_log(:,1),'omitnan');
APA_bacterial_log_std=std(APA_bacterial_log(:,1),'omitnan');
APA_phytoplankton_log_std=std(APA_phytoplankton_log(:,1),'omitnan');
APA_trichodesmium_log_std=std(APA_trichodesmium_log(:,1),'omitnan');
APA_log_std=[APA_bulk_log_std,APA_dissolved_log_std,APA_particulate_log_std,APA_bacterial_log_std,APA_phytoplankton_log_std,APA_trichodesmium_log_std];
% calculate number of data points
APA_bulk_log_num=length(find(~isnan(APA_bulk_log)));
APA_dissolved_log_num=length(find(~isnan(APA_dissolved_log)));
APA_particulate_log_num=length(find(~isnan(APA_particulate_log)));
APA_bacterial_log_num=length(find(~isnan(APA_bacterial_log)));
APA_phytoplankton_log_num=length(find(~isnan(APA_phytoplankton_log)));
APA_trichodesmium_log_num=length(find(~isnan(APA_trichodesmium_log)));
APA_log_number=[APA_bulk_log_num,APA_dissolved_log_num,APA_particulate_log_num,APA_bacterial_log_num,APA_phytoplankton_log_num,APA_trichodesmium_log_num];
%%
% calculate the critical value
for i=1:length(APA_log_number)
    P(i)=1-1./(4.*APA_log_number(i));%#
end 
APA_log_criterion=norminv(P,APA_log_mean,APA_log_std);% APA_log_criterion is the critical value to identify outliers
%%
% MUF-P substrate
% APA bulk fraction
% 'n' variable will save APA data after flagged outliers
clear i
i=1;
n_bulk = nan(length(APA_bulk_log),1); % n_bulk contains real APA values that are not identified as outliers, and is 0 for an identified outlier.
for t=1:length(APA_bulk_log)
    if APA_bulk_log(t)>APA_log_criterion(1)
        n_bulk(t)=0;
  		num_bulk(i)=APA_bulk_log(t);% num_bulk is the outlier value when an outlier is identified
  		i=i+1;
	else
  		n_bulk(t)=APA_bulk_log(t);
    end
end
k_bulk=sum(n_bulk==0);% number of outliers 
n_m_bulk=nanmean(n_bulk(n_bulk~=0)); % calculte the mean of log10-transform APA values after removing the identified outliers.
n_std_bulk=nanstd(n_bulk(n_bulk~=0)); % calculte standard deviation of log10-transform APA values after removing the identified outliers.
%%
% APA dissolved fraction
clear i t
i=1;
n_dissolved = nan(length(APA_dissolved_log),1); % n_dissolved contains real log10-transformed APA values that are not identified as outliers, and is 0 for an identified outlier.
for t=1:length(APA_dissolved_log)
	if APA_dissolved_log(t)>APA_log_criterion(2)
  		n_dissolved(t)=0;
        num_dissolved(i)=APA_dissolved_log(t);%num_dissolved is the outlier value when an outlier is identified
  		i=i+1;
	else
  		n_dissolved(t)=APA_dissolved_log(t);
  	end
end
k_dissolved=sum(n_dissolved==0);% number of outliers
n_m_dissolved=nanmean(n_dissolved(n_dissolved~=0)); %calculte the mean of log10-transform APA values after removing the identified outliers.
n_std_dissolved=nanstd(n_dissolved(n_dissolved~=0));%calculte standard deviation of log10-transform APA values after removing the identified outliers
%%
% APA particulate fraction
clear i t
i=1;
n_particulate = nan(length(APA_particulate_log),1); % n_particulate contains real log10-transformed APA values that are not identified as outliers, and is 0 for an identified outlier.
for t=1:length(APA_particulate_log)
	if APA_particulate_log(t)>APA_log_criterion(3)%% 
  		n_particulate(t)=0;
        num_particulate(i)=APA_particulate_log(t);%num_particulate is the outlier value when an outlier is identified
  		i=i+1;
	else
  		n_particulate(t)=APA_particulate_log(t);
  	end
end
k_particulate=sum(n_particulate==0);% number of outliers
n_m_particulate=nanmean(n_particulate(n_particulate~=0)); %calculte the mean of log10-transform APA values after removing the identified outliers.
n_std_particulate=nanstd(n_dissolved(n_particulate~=0));%calculte standard deviation of log10-transform APA values after removing the identified outliers
%%

%%
% APA bacterial fraction
clear i t
i=1;
n_bacterial = nan(length(APA_bacterial_log),1); % n_bacterial contains real log10-transformed APA values that are not identified as outliers, and is 1 for an identified outlier.
for t=1:length(APA_bacterial_log)
	if APA_bacterial_log(t)>APA_log_criterion(4)
  		n_bacterial(t)=1;%%n_bacterial is 1 for an identified outlier. 
        num_bacterial(i)=APA_bacterial_log(t);%num_bacterial is the outlier value when an outlier is identified
  		i=i+1;
	else
  		n_bacterial(t)=APA_bacterial_log(t);
  	end
end
k_bacterial=sum(n_bacterial==1);% number of outliers
n_m_bacterial=nanmean(n_bacterial(n_bacterial~=1)); %calculte the mean of log10-transform APA values after removing the identified outliers.
n_std_bacterial=nanstd(n_bacterial(n_bacterial~=1));%calculte standard deviation of log10-transform APA values after removing the identified outliers

%%
% APA phytoplankton fraction
clear i t
i=1;
n_phytoplankton = nan(length(APA_phytoplankton_log),1); % n_phytoplankton contains real log10-transformed APA values that are not identified as outliers, and is 0 for an identified outlier.
for t=1:length(APA_phytoplankton_log)
	if APA_phytoplankton_log(t)>APA_log_criterion(5)
  		n_phytoplankton(t)=0;%%n_phytoplankton is 0 for an identified outlier. 
        num_phytoplankton(i)=APA_phytoplankton_log(t);%num_phytoplankton is the outlier value when an outlier is identified
  		i=i+1;
	else
  		n_phytoplankton(t)=APA_phytoplankton_log(t);
  	end
end
k_phytoplankton=sum(n_phytoplankton==0);% number of outliers
n_m_phytoplankton=nanmean(n_phytoplankton(n_phytoplankton~=0)); %calculte the mean of log10-transform APA values after removing the identified outliers.
n_std_phytoplankton=nanstd(n_phytoplankton(n_phytoplankton~=0));%calculte standard deviation of log10-transform APA values after removing the identified outliers

%%% APA trichodesmium fraction
clear i t
i=1;
n_trichodesmium = nan(length(APA_trichodesmium_log),1); % n_trichodesmium contains real log10-transformed APA values that are not identified as outliers, and is 0 for an identified outlier.
for t=1:length(APA_trichodesmium_log)
	if APA_trichodesmium_log(t)>APA_log_criterion(6)
  		n_trichodesmium(t)=0;%%n_trichodesmium is 0 for an identified outlier. 
        num_trichodesmium(i)=APA_trichodesmium_log(t);%num_trichodesmium is the outlier value when an outlier is identified
  		i=i+1;
	else
  		n_trichodesmium(t)=APA_trichodesmium_log(t);
  	end
end
k_trichodesmium=sum(n_trichodesmium==0);% number of outliers
n_m_trichodesmium=nanmean(n_trichodesmium(n_trichodesmium~=0)); %calculte the mean of log10-transform APA values after removing the identified outliers.
n_std_trichodesmium=nanstd(n_trichodesmium(n_trichodesmium~=0));%calculte standard deviation of log10-transform APA values after removing the identified outliers


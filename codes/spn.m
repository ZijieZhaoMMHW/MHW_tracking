%% radius
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% radius_sp=cell(length(tracks),1);
% [lat2,lon2]=meshgrid(lat_used,lon_used);
% 
% t_full=tracks_new(1,:);
% tmin=nanmin(cellfun(@nanmin,t_full));
% tmax=nanmax(cellfun(@nanmax,t_full));
% 
% for i=1:length(tracks)
%     tic
%     t_here=tracks_new{1,i};
%     radius_s=cell(length(t_here),1);
%     xloc=tracks_new{2,i};
%     yloc=tracks_new{3,i};
%     for j=1:length(t_here)
% %         mhw_here=sst_anom(:,:,t_here(j));
%         ind_here=sub2ind([360 180],xloc{j},yloc{j});
%         img_here=zeros(360,180);
%         img_here(ind_here)=1;
%         bw=bwconncomp(img_here,8);
%         %sst_b=NaN(100,100,bw.NumObjects);
%         radius_b=NaN(1,1);
%         for b=1:1
%             idx_here=ind_here;
%             center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
% %            center_b(b,:)=center_here;
%             dist_here=geodist(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%             
% %             area_here=length(idx_here)*0.0625;
% %             area_b(b)=area_here;
%             
%             radius_b(b)=nanmax(dist_here);
%         end
% %         sst_s(:,:,j)=nansum(sst_b.*reshape(area_b,1,1,length(area_b)),3)./nansum(area_b);
%         %center_s(j,:)=nansum(center_b.*area_b,1)./nansum(area_b);
% %         x_s{j}=x_b;y_s{j}=y_b;
%         %area_s{j}=area_b;
%         radius_s{j}=radius_b;
%     end
% %     sst_sp{i}=sst_s;
% %     x_sp{i}=x_s;
% %     y_sp{i}=y_s;
% %     area_sp{i}=area_s;
%     radius_sp{i}=radius_s;
%     toc
% end
% radius_s=NaN(length(radius_sp),1);
% for i=1:length(radius_sp)
%     radius_here=radius_sp{i};
%     radius_here=cellfun(@nanmax,radius_here);
%     radius_s(i)=nanmax(radius_here);
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)],'radius_s','-v7.3'); 

%% spn
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'temp','q','sp','tcc','tp','ws'};
% data_anom=NaN(360,180,length(t_min:t_max),6);
% 
% for i=1:6
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),6);
% 
% for m=1:6
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         tic
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%         toc
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_all' num2str(d)],'sst_spn_full','-v7.3');

%% spn-z500
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'z500'};
% data_anom=NaN(360,180,length(t_min:t_max),1);
% 
% for i=1:1
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),1);
% 
% for m=1:1
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         tic
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%         toc
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_z500_' num2str(d)],'sst_spn_full','-v7.3');

%% spn-2 radius
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'temp','q','sp','tcc','tp','z500'};
% data_anom=NaN(360,180,length(t_min:t_max),6);
% 
% for i=1:6
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),6);
% 
% for m=1:6
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         tic
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%         toc
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_all_2r' num2str(d)],'sst_spn_full','-v7.3');


%% spn-2 radius-all
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'temp','q','sp','tcc','tp','z500'};
% data_anom=NaN(360,180,length(t_min:t_max),6);
% 
% for i=1:6
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),6);
% 
% for m=1:6
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         tic
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%         toc
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_all_2r_raw' num2str(d)],'sst_spn_full','-v7.3');

%% spn-2 radius-all-FX
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'FXp','FXt','FYp','FYt','FXpt','FYpt'};
% data_anom=NaN(360,180,length(t_min:t_max),6);
% 
% for i=1:6
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),6);
% 
% for m=1:6
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         tic
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%         toc
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_f_2r_raw' num2str(d)],'sst_spn_full','-v7.3');

%% spn-2 radius-all-FX-anom
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'q','latent','sensible','short','long'};
% data_anom=NaN(360,180,length(t_min:t_max),5);
% 
% for i=1:5
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),5);
% 
% for m=1:5
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_f_2r_q_dhw' num2str(d)],'sst_spn_full','-v7.3');

%% direction
% load('hwmax_metricts');
% 
% directions=NaN(length(area_ts),5,2);
% 
% for i=1:length(area_ts)
%     tracks_here=tracks_line{i};
%     t_here=(1.5:(size(tracks_here,1)-0.5))';
%     d_here=NaN(size(tracks_here,1)-1,2);
%     for j=1:(size(tracks_here,1)-1)
%         x2 = tracks_here(j,1);
%         y2 = tracks_here(j,2);
%         
%         x1 = tracks_here(j+1,1);
%         y1 = tracks_here(j+1,2);
%         
%         [~,x,y]=lldistkm(y1,y2,x1,x2);
%         
%         PC1=x;
%         PC2=y;
%         
% %         d_here(j)=(angle(PC2(:)+1i*PC1(:)))./pi*180;
%         d_here(j,:)=[PC1 PC2];
%     end
%     
% end

%% 
% d= depth_s;
% addpath /g/data/v45/zz6006/mhw_timescale/CDT-master/cdt
% addpath /g/data/v45/zz6006/mhw_timescale/CDT-master/cdt/cdt_data
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% load('hwmax_metricts');
% tracks_line=tracks_line(idx_here);
% dist2co=NaN(length(tracks_line),1);
% 
% for i=1:length(tracks_line)
%     tracks_here=tracks_line{i};
%     o_here=tracks_here(end,:);
%     dist2co(i)=(-1).^double(island(o_here(2),o_here(1))).*dist2coast(o_here(2),o_here(1));
% end
% 
% for i=1:length(tracks_line)
%     tracks_here=tracks_line{i};
%     o_here=tracks_here(end,:);
%     lats(i)=o_here(2);
% end
% save(['dis2coe' num2str(d)],'dist2co','lats');

%% spn-2 radius-all-w850 full/anomaly
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'w850','w850_anom'};
% data_anom=NaN(360,180,length(t_min:t_max),length(name_used));
% 
% for i=1:length(name_used)
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     if i==1
%     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
%     else
%         eval(['anom_here=mfile.data_anom(:,:,t_min:t_max);'])
%     end
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),length(name_used));
% 
% for m=1:length(name_used)
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_2r_w850_' num2str(d)],'sst_spn_full','-v7.3');

%% spn-2 radius-all-w500 full/anomaly
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'w500','w500_anom'};
% data_anom=NaN(360,180,length(t_min:t_max),length(name_used));
% 
% for i=1:length(name_used)
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     if i==1
%     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
%     else
%         eval(['anom_here=mfile.data_anom(:,:,t_min:t_max);'])
%     end
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),length(name_used));
% 
% for m=1:length(name_used)
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_2r_w500_' num2str(d)],'sst_spn_full','-v7.3');


%% spn-2 radius-all-FX-anom
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'q','latent','sensible','short','long'};
% data_anom=NaN(360,180,length(t_min:t_max),5);
% 
% for i=1:5
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),5);
% 
% for m=1:5
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         sst_here=sst_sp{i};
%         tb=(1:size(sst_here,3))./size(sst_here,3);
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_2r_q_raw' num2str(d)],'sst_spn_full','-v7.3');

%% spn ts
% d= depth_s;
% lona=0.5:1:360;
% lata=-89.5:1:89.5;
% [lata,lona]=meshgrid(lata,lona);
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'sp','sensible','w850','wy'};
% data_anom=NaN(360,180,length(t_min:t_max),length(name_used));
% load('/g/data/v45/zz6006/mmhw_cmip/data_re/used_loc')
% for i=1:length(name_used)
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     if i<=2
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     else
%         eval(['anom_here=mfile.data_anom(:,:,t_min:t_max);'])
%     end
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% % name_used={'sp'};
% % 
% % for i=1:1
% %     tic
% %     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
% %     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
% %     data_anom(:,:,:,10)=anom_here;
% %     toc
% % end
% load('/scratch/v45/zz6006/mmhw_cmip/area_used')
% q_ts=cell(length(tracks),length(name_used));
% for m=1:4
%     q_anom=data_anom(:,:,:,m);
%     for i=1:length(tracks)
%         tic
%         t_here=tracks_new{1,i};
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         q_here=NaN(length(t_here),1);
%         for j=1:length(t_here)
%             mhw_here=q_anom(:,:,t_here(j)-t_min+1);
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             q_here(j)=nansum(mhw_here(ind_here).*area_used(ind_here));
%         end
%         q_ts{i,m}=q_here;
%         toc
%     end
% end
% 
% % for i=1:length(tracks)
% %     tic
% %     t_here=tracks_new{1,i};
% %     xloc=tracks_new{2,i};
% %     yloc=tracks_new{3,i};
% %     q_here=NaN(length(t_here),1);
% %     for j=1:length(t_here)
% %         ind_here=sub2ind([360 180],xloc{j},yloc{j});
% %         q_here(j)=nansum(cosd(lata(ind_here)).*land_idx(ind_here))./...
% %         nansum(cosd(lata(ind_here)));
% %     end
% %     q_ts{i,11}=q_here;
% %     toc
% % end
% 
% save(['q_ts' num2str(d)],'q_ts','-v7.3');

%% q_ts_norm
% for d=1:200
%     tic
%     load(['q_ts' num2str(d)]);
%     q_ts_norm=NaN(size(q_ts,1),4,25);
%     for i=1:size(q_ts,1)
%         for j=1:size(q_ts,2)
%             tsb=q_ts{i,j};
%             tsb=tsb./nanmax(abs(tsb));
%             tb=linspace(0,1,length(tsb));
%             ta=linspace(0,1,25);
%             tsa=interp1(tb,tsb,ta,'linear','extrap');
%             q_ts_norm(i,j,:)=tsa;
%         end
%     end
%     save(['q_ts_norm' num2str(d)],'q_ts_norm','-v7.3');
%     toc
% end

%% spn-2 radius-all
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'sp'};
% data_anom=NaN(360,180,length(t_min:t_max),length(name_used));
% 
% for i=1:length(name_used)
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),length(name_used));
% 
% for m=1:length(name_used)
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/sp_sp' num2str(d)],'sst_sp','-v7.3');

%%
% d= depth_s;
% load(['/g/data/v45/zz6006/mhw_hw_mhw/sp_sp' num2str(d)]);
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% sp_ts=cell(100,1);
% res=50;
% x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
% radius_range=linspace(1,0,res);
% radius_range=radius_range(end:-1:1);
% for r=1:res
%     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
% end
% x_full=x_full(:);
% y_full=y_full(:);
% for i=1:length(sst_sp)
%     sp_here=sst_sp{i};
%     sp_here=reshape(sp_here,50*50,size(sp_here,3));
%     inner_here=sp_here(sqrt(x_full.^2+y_full.^2)<=0.5,:);
%     outer_here=sp_here(sqrt(x_full.^2+y_full.^2)>0.5,:);
%     sp_ts{i}=nanmean(inner_here)-nanmean(outer_here);
% end
% save(['sp_ts' num2str(d)],'sp_ts');

%% combined metric
% d= depth_s;
% lona=0.5:1:360;
% lata=-89.5:1:89.5;
% [lata,lona]=meshgrid(lata,lona);
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'sp','sensible','w850','wy'};
% data_anom=NaN(360,180,length(t_min:t_max),length(name_used));
% load('/g/data/v45/zz6006/mmhw_cmip/data_re/used_loc')
% for i=1:length(name_used)
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     if i<=2
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     else
%         eval(['anom_here=mfile.data_anom(:,:,t_min:t_max);'])
%     end
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% % name_used={'sp'};
% % 
% % for i=1:1
% %     tic
% %     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
% %     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
% %     data_anom(:,:,:,10)=anom_here;
% %     toc
% % end
% load('/scratch/v45/zz6006/mmhw_cmip/area_used')
% q_ts=cell(length(tracks),length(name_used));
% for m=1:4
%     q_anom=data_anom(:,:,:,m);
%     for i=1:length(tracks)
%         tic
%         t_here=tracks_new{1,i};
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         q_here=NaN(length(t_here),1);
%         for j=1:length(t_here)
%             mhw_here=q_anom(:,:,t_here(j)-t_min+1);
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             q_here(j)=nansum(mhw_here(ind_here).*area_used(ind_here));
%         end
%         q_ts{i,m}=q_here;
%         toc
%     end
% end
% load(['sp_ts' num2str(d)]);
% q_ts=cat(2,q_ts(:,2:end),sp_ts);
% 
% q_ts_norm=NaN(size(q_ts,1),4);
% for i=1:size(q_ts,1)
%     for j=1:size(q_ts,2)
%         tsb=q_ts{i,j};
%         tb=linspace(0,1,length(tsb));
%         ta=linspace(0,1,25);
%         ta=ta(ta>=0.25 & ta<=0.75);
%         tsa=interp1(tb,tsb,ta,'linear','extrap');
%         q_ts_norm(i,j)=nanmean(tsa);
%     end
% end
% save(['subjts' num2str(d)],'q_ts_norm','-v7.3');

%% cluster spn
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'temp','sensible','sp','w850','wy','sp'};
% data_anom=NaN(360,180,length(t_min:t_max),5);
% 
% for i=1:length(name_used)
%     tic
%     if i<4
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     elseif i<=5
%         mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     eval(['anom_here=mfile.data_anom(:,:,t_min:t_max);'])
%     else
%         mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
%     end
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),length(name_used));
% 
% for m=1:length(name_used)
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         sst_here=sst_sp{i};
%         %tb=(1:size(sst_here,3))./size(sst_here,3);
%         tb=linspace(0,1,size(sst_here,3));
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_sub' num2str(d)],'sst_spn_full','-v7.3');

%% ts hw
% d= depth_s;
% lona=0.5:1:360;
% lata=-89.5:1:89.5;
% [lata,lona]=meshgrid(lata,lona);
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'sp','sensible','w850','wy'};
% data_anom=NaN(360,180,length(t_min:t_max),length(name_used));
% load('/g/data/v45/zz6006/mmhw_cmip/data_re/used_loc')
% for i=1:length(name_used)
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     if i<=2
%     eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
%     else
%         eval(['anom_here=mfile.data_anom(:,:,t_min:t_max);'])
%     end
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% % name_used={'sp'};
% % 
% % for i=1:1
% %     tic
% %     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
% %     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
% %     data_anom(:,:,:,10)=anom_here;
% %     toc
% % end
% load('/scratch/v45/zz6006/mmhw_cmip/area_used')
% q_ts=cell(length(tracks),length(name_used));
% for m=1:4
%     q_anom=data_anom(:,:,:,m);
%     for i=1:length(tracks)
%         tic
%         t_here=tracks_new{1,i};
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         q_here=NaN(length(t_here),1);
%         for j=1:length(t_here)
%             mhw_here=q_anom(:,:,t_here(j)-t_min+1);
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             q_here(j)=nansum(mhw_here(ind_here).*area_used(ind_here));
%         end
%         q_ts{i,m}=q_here;
%         toc
%     end
% end
% load(['sp_ts' num2str(d)]);
% q_ts=cat(2,q_ts(:,2:end),sp_ts);
% 
% q_ts_norm=NaN(size(q_ts,1),25,4);
% for i=1:size(q_ts,1)
%     for j=1:size(q_ts,2)
%         tsb=q_ts{i,j};
%         tb=linspace(0,1,length(tsb));
%         ta=linspace(0,1,25);
%         %ta=ta(ta>=0.25 & ta<=0.75);
%         tsa=interp1(tb,tsb,ta,'linear','extrap');
%         q_ts_norm(i,:,j)=tsa;
%     end
% end
% save(['subjts_all' num2str(d)],'q_ts_norm','-v7.3');


%% ts hw m
% d= depth_s;
% lona=0.5:1:360;
% lata=-89.5:1:89.5;
% [lata,lona]=meshgrid(lata,lona);
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'wx'};
% data_anom=NaN(360,180,length(t_min:t_max),length(name_used));
% load('/g/data/v45/zz6006/mmhw_cmip/data_re/used_loc')
% for i=1:length(name_used)
%     tic
%     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
%     
%         eval(['anom_here=mfile.wx_full(:,:,t_min:t_max);'])
%   
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% % name_used={'sp'};
% % 
% % for i=1:1
% %     tic
% %     mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_full']);
% %     eval(['anom_here=mfile.' name_used{i} '_full(:,:,t_min:t_max);'])
% %     data_anom(:,:,:,10)=anom_here;
% %     toc
% % end
% load('/scratch/v45/zz6006/mmhw_cmip/area_used')
% q_ts=cell(length(tracks),length(name_used));
% data_anom=data_anom-nanmean(data_anom);
% for m=1:1
%     q_anom=data_anom(:,:,:,m);
%     for i=1:length(tracks)
%         tic
%         t_here=tracks_new{1,i};
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         q_here=NaN(length(t_here),1);
%         for j=1:length(t_here)
%             mhw_here=q_anom(:,:,t_here(j)-t_min+1);
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             q_here(j)=nansum(mhw_here(ind_here).*area_used(ind_here))./nansum(area_used(ind_here));
%         end
%         q_ts{i,m}=q_here;
%         toc
%     end
% end
% 
% q_ts_norm=NaN(size(q_ts,1),25);
% for i=1:size(q_ts,1)
%     for j=1:size(q_ts,2)
%         tsb=q_ts{i,j};
%         tb=linspace(0,1,length(tsb));
%         ta=linspace(0,1,25);
%         %ta=ta(ta>=0.25 & ta<=0.75);
%         tsa=interp1(tb,tsb,ta,'linear','extrap');
%         q_ts_norm(i,:,j)=tsa;
%     end
% end
% save(['subjts_wx' num2str(d)],'q_ts_norm','-v7.3');

%% cluster spn
% d= depth_s;
% if d<200
%     idx_here=(1:100)+(d-1)*100;
% else
%     idx_here=19901:20023;
% end
% 
% load(['/scratch/v45/zz6006/mmhw_cmip/area_used']);
% mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
% load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
% tracks=mfile.tracks(:,idx_here);
% 
% tracks_new=squeeze(struct2cell(tracks));
% t_used=tracks_new(1,:);
% t_min=nanmin(cellfun(@nanmin,t_used));
% t_max=nanmax(cellfun(@nanmax,t_used));
% 
% name_used={'wx'};
% data_anom=NaN(360,180,length(t_min:t_max),5);
% 
% for i=1:length(name_used)
%     tic
%     
%         mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
%     eval(['anom_here=mfile.data_anom(:,:,t_min:t_max);'])
%     
%     data_anom(:,:,:,i)=anom_here;
%     toc
% end
% 
% lon_used=0.5:1:360;
% lat_used=-89.5:89.5;
% 
% %addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
% sst_spn_full=NaN(50,50,5,length(tracks),length(name_used));
% 
% for m=1:length(name_used)
%     tic
%     sst_anom=data_anom(:,:,:,m);
%     sst_sp=cell(length(tracks),1);
%     x_sp=cell(length(tracks),1);
%     y_sp=cell(length(tracks),1);
%     area_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     res=50;
%     
%     for i=1:length(tracks)
%         t_here=tracks_new{1,i}-t_min+1;
%         sst_s=NaN(res,res,length(t_here));
%         x_s=cell(length(t_here),1);
%         y_s=cell(length(t_here),1);
%         center_s=NaN(length(t_here),2);
%         area_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([360 180],xloc{j},yloc{j});
%             img_here=zeros(360,180);
%             img_here(ind_here)=1;
%             %bw=bwconncomp(img_here,8);
%             sst_b=NaN(res,res,1);
%             x_b=NaN(res,res,1);
%             y_b=NaN(res,res,1);
%             %         area_b=NaN(1,1);
%             center_b=NaN(1,2);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 center_b(b,:)=center_here;
%                 %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
%                 %             +(lat2(idx_here)-center_here(2)).^2);
%                 %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 radius_here=2*radius_s(i);
%                 radius_range=linspace(radius_here,0,res);
%                 radius_range=radius_range(end:-1:1);
%                 x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
%                 for r=1:res
%                     x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
%                     y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
%                 end
%                 [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
%                 F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
%                 sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
%                 k=boundary(lond(idx_here),latd(idx_here),1);
%                 if isempty(k)
%                     k=1:length(lond(idx_here));
%                 end
%                 lon2k=lond(idx_here);lat2k=latd(idx_here);
%                 [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
%                 %sst_b_here(~in & ~on)=0;
%                 sst_b(:,:,b)=sst_b_here;
%                 x_b(:,:,b)=x_full;
%                 y_b(:,:,b)=y_full;
%             end
%             sst_s(:,:,j)=sst_b;
%             center_s(j,:)=center_b;
%             x_s{j}=x_b;y_s{j}=y_b;
%         end
%         sst_sp{i}=sst_s;
%         x_sp{i}=x_s;
%         y_sp{i}=y_s;
%     end
%     
%     t_used=0:0.25:1;
%     radius_range=linspace(1,0,50);
%     radius_range=radius_range(end:-1:1);
%     x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
%     for r=1:50
%         x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
%         y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
%     end
%     
%     sst_spn=NaN(50,50,5,length(sst_sp));
%     ta=(0:0.25:1)';
%     for i=1:length(sst_sp)
%         sst_here=sst_sp{i};
%         %tb=(1:size(sst_here,3))./size(sst_here,3);
%         tb=linspace(0,1,size(sst_here,3));
%         sst_b=NaN(50,50,length(ta));
%         for x=1:size(sst_here,1)
%             for y=1:size(sst_here,2)
%                 ts_b=squeeze(sst_here(x,y,:));
%                 ts_a=interp1(tb,ts_b,ta,'linear','extrap');
%                 sst_b(x,y,:)=ts_a;
%             end
%         end
%         sst_spn(:,:,:,i)=sst_b;
%     end
%     sst_spn_full(:,:,:,:,m)=sst_spn;
%     
%     toc
% end
% save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_wx' num2str(d)],'sst_spn_full','-v7.3');

%% cluster spn
li=15877;
lon_used=0.5:1:360;
lat_used=-89.5:89.5;
load('tracks_draw');
load(['radius_f']);
radius_s=radius_f(li);

tracks_new=squeeze(struct2cell(tracks));
t_used=tracks_new(1,:);
t_min=nanmin(cellfun(@nanmin,t_used));
t_max=nanmax(cellfun(@nanmax,t_used));

load('temp_sub');
data_anom=temp_sub;

lon_used=0.5:1:360;
lat_used=-89.5:89.5;

%addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
sst_spn_full=NaN(100,100,5,length(tracks));

for m=1:1
    tic
    sst_anom=data_anom(:,:,:,m);
    sst_sp=cell(length(tracks),1);
    x_sp=cell(length(tracks),1);
    y_sp=cell(length(tracks),1);
    area_sp=cell(length(tracks),1);
    [lat2,lon2]=meshgrid(lat_used,lon_used);
    res=100;
    
    for i=1:length(tracks)
        t_here=tracks_new{1,i}-t_min+1;
        sst_s=NaN(res,res,length(t_here));
        x_s=cell(length(t_here),1);
        y_s=cell(length(t_here),1);
        center_s=NaN(length(t_here),2);
        area_s=cell(length(t_here),1);
        xloc=tracks_new{2,i};
        yloc=tracks_new{3,i};
        for j=[3 436 530];
            mhw_here=sst_anom(:,:,t_here(j));
            ind_here=sub2ind([360 180],xloc{j},yloc{j});
            img_here=zeros(360,180);
            img_here(ind_here)=1;
            %bw=bwconncomp(img_here,8);
            sst_b=NaN(res,res,1);
            x_b=NaN(res,res,1);
            y_b=NaN(res,res,1);
            %         area_b=NaN(1,1);
            center_b=NaN(1,2);
            for b=1:1
                idx_here=ind_here;
                center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
                center_b(b,:)=center_here;
                %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
                %             +(lat2(idx_here)-center_here(2)).^2);
                %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
                [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
                %             area_here=length(idx_here)*0.0625;
                %             area_b(b)=area_here;
                radius_here=2*radius_s(i)./3.5;
                radius_range=linspace(radius_here,0,res);
                radius_range=radius_range(end:-1:1);
                x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
                for r=1:res
                    x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
                    y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
                end
                [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
                F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
                sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
                k=boundary(lond(idx_here),latd(idx_here),1);
                if isempty(k)
                    k=1:length(lond(idx_here));
                end
                lon2k=lond(idx_here);lat2k=latd(idx_here);
                [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
                sst_b_here(~in & ~on)=0;
                sst_b(:,:,b)=sst_b_here;
                x_b(:,:,b)=x_full;
                y_b(:,:,b)=y_full;
            end
            sst_s(:,:,j)=sst_b;
            center_s(j,:)=center_b;
            x_s{j}=x_b;y_s{j}=y_b;
        end
        sst_sp{i}=sst_s;
        x_sp{i}=x_s;
        y_sp{i}=y_s;
    end
    
    t_used=0:0.25:1;
    radius_range=linspace(1,0,50);
    radius_range=radius_range(end:-1:1);
    x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
    for r=1:50
        x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
        y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
    end
    
    sst_spn=NaN(50,50,5,length(sst_sp));
    
    
    toc
end
save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_tpsh_r' num2str(d)],'sst_spn_full','-v7.3');

%%
d= depth_s;
if d<200
    idx_here=(1:100)+(d-1)*100;
else
    idx_here=19901:20023;
end

load('area_used');
mfile=matfile('/scratch/v45/zz6006/mhw_hw_mhw/hw/tracks_hw_max');
load(['/g/data/v45/zz6006/mhw_hw_mhw/radius_max' num2str(d)]);
tracks=mfile.tracks(:,idx_here);

tracks_new=squeeze(struct2cell(tracks));
t_used=tracks_new(1,:);
t_min=nanmin(cellfun(@nanmin,t_used));
t_max=nanmax(cellfun(@nanmax,t_used));

name_used={'temp','q','sp','tcc','tp','ws'};
data_anom=NaN(360,180,length(t_min:t_max),6);

for i=1:6
    tic
    mfile=matfile(['/g/data/v45/zz6006/mhw_hw_mhw/' name_used{i} '_anom_full']);
    eval(['anom_here=mfile.' name_used{i} '_anom_full(:,:,t_min:t_max);'])
    data_anom(:,:,:,i)=anom_here;
    toc
end

lon_used=0.5:1:360;
lat_used=-89.5:89.5;

addpath /scratch/v45/zz6006/mhw_dynamical/mhw_new
sst_spn_full=NaN(50,50,5,length(tracks),6);

for m=1:6
    tic
    sst_anom=data_anom(:,:,:,m);
    sst_sp=cell(length(tracks),1);
    x_sp=cell(length(tracks),1);
    y_sp=cell(length(tracks),1);
    area_sp=cell(length(tracks),1);
    [lat2,lon2]=meshgrid(lat_used,lon_used);
    res=50;
    
    for i=1:length(tracks)
        t_here=tracks_new{1,i}-t_min+1;
        sst_s=NaN(res,res,length(t_here));
        x_s=cell(length(t_here),1);
        y_s=cell(length(t_here),1);
        center_s=NaN(length(t_here),2);
        area_s=cell(length(t_here),1);
        xloc=tracks_new{2,i};
        yloc=tracks_new{3,i};
        for j=1:length(t_here)
            mhw_here=sst_anom(:,:,t_here(j));
            ind_here=sub2ind([360 180],xloc{j},yloc{j});
            img_here=zeros(360,180);
            img_here(ind_here)=1;
            %bw=bwconncomp(img_here,8);
            sst_b=NaN(res,res,1);
            x_b=NaN(res,res,1);
            y_b=NaN(res,res,1);
            %         area_b=NaN(1,1);
            center_b=NaN(1,2);
            for b=1:1
                idx_here=ind_here;
                center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
                center_b(b,:)=center_here;
                %             dist_here=sqrt((lon2(idx_here)-center_here(1)).^2 ...
                %             +(lat2(idx_here)-center_here(2)).^2);
                %            dist_here=lldistkm(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
                [~,lond,latd]=geodist(center_here(2),lat2,center_here(1),lon2);
                %             area_here=length(idx_here)*0.0625;
                %             area_b(b)=area_here;
                radius_here=radius_s(i);
                radius_range=linspace(radius_here,0,res);
                radius_range=radius_range(end:-1:1);
                x_full=NaN(res,res);y_full=NaN(res,res);%radius -data
                for r=1:res
                    x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,res));
                    y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,res));
                end
                [in,on]=inpolygon(lond(:),latd(:),x_full(end,:),y_full(end,:));
                F=scatteredInterpolant(lond(in | on),latd(in | on),mhw_here(in | on),'linear','linear');
                sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
                k=boundary(lond(idx_here),latd(idx_here),1);
                if isempty(k)
                    k=1:length(lond(idx_here));
                end
                lon2k=lond(idx_here);lat2k=latd(idx_here);
                [in,on]=inpolygon(x_full(:),y_full(:),lon2k(k),lat2k(k));
                %sst_b_here(~in & ~on)=0;
                sst_b(:,:,b)=sst_b_here;
                x_b(:,:,b)=x_full;
                y_b(:,:,b)=y_full;
            end
            sst_s(:,:,j)=sst_b;
            center_s(j,:)=center_b;
            x_s{j}=x_b;y_s{j}=y_b;
        end
        sst_sp{i}=sst_s;
        x_sp{i}=x_s;
        y_sp{i}=y_s;
    end
    
    t_used=0:0.25:1;
    radius_range=linspace(1,0,50);
    radius_range=radius_range(end:-1:1);
    x_full=NaN(50,50);y_full=NaN(50,50);%radius -data
    for r=1:50
        x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,50));
        y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,50));
    end
    
    sst_spn=NaN(50,50,5,length(sst_sp));
    ta=(0:0.25:1)';
    for i=1:length(sst_sp)
        tic
        sst_here=sst_sp{i};
        tb=(1:size(sst_here,3))./size(sst_here,3);
        sst_b=NaN(50,50,length(ta));
        for x=1:size(sst_here,1)
            for y=1:size(sst_here,2)
                ts_b=squeeze(sst_here(x,y,:));
                ts_a=interp1(tb,ts_b,ta,'linear','extrap');
                sst_b(x,y,:)=ts_a;
            end
        end
        sst_spn(:,:,:,i)=sst_b;
        toc
    end
    sst_spn_full(:,:,:,:,m)=sst_spn;
    
    toc
end
save(['/g/data/v45/zz6006/mhw_hw_mhw/spn_full_all' num2str(d)],'sst_spn_full','-v7.3');


%%
idx_used=[3 436 530];

radius_range=linspace(1,0,100);
radius_range=radius_range(end:-1:1);
x_full=NaN(100,100);y_full=NaN(100,100);%radius -data
for r=1:100
    x_full(r,:)=radius_range(r)*sin(linspace(0,2*pi,100));
    y_full(r,:)=radius_range(r)*cos(linspace(0,2*pi,100));
end

for i=1:3
    subplot(3,1,i);
    pcolor(x_full,y_full,sst_sp{1}(:,:,idx_used(i)));
    colormap(cmocean('amp'));
    hold on
    plot(x_full(end,:),y_full(end,:),'k','linewidth',2);
    hold on
    plot(x_full(50,:),y_full(50,:),'k','linewidth',2);
    shading interp
    hold on
    caxis([0 3]);
    axis equal
        axis tight
        box on
        axis off
        set(gca,'linewidth',2);
end
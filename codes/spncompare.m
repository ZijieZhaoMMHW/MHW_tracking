%% radius - sc
% d = depth_s;
% if d<191
% idx_here=(1:5)+(d-1)*5;
% else
%     idx_here=951:958;
% end
% 
% load('areas');
% area_used=areas;
% mfile=matfile('tracks3d5r');
% tracks=mfile.tracks(:,idx_here);
% load('lonlat');
% tracks_new=squeeze(struct2cell(tracks));
% lon_used=xh;
% lat_used=yh;
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
% %         hw_here=sst_anom(:,:,t_here(j));
%         ind_here=sub2ind([775 845],xloc{j},yloc{j});
%         img_here=zeros(775,845);
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
% save(['radius_sc' num2str(d)],'radius_s','-v7.3');

%% spn temp
d = depth_s;
if d<191
idx_here=(1:5)+(d-1)*5;
else
    idx_here=951:958;
end

load('areas');
area_used=areas;
mfile=matfile('tracks3d5r');
tracks=mfile.tracks(:,idx_here);
load('lonlat');
tracks_new=squeeze(struct2cell(tracks));
lon_used=xh;
lat_used=yh;

sst_sp=cell(length(tracks),1);
x_sp=cell(length(tracks),1);
y_sp=cell(length(tracks),1);
area_sp=cell(length(tracks),1);
[lat2,lon2]=meshgrid(lat_used,lon_used);
res=50;

t_full=tracks_new(1,:);
t_min=nanmin(cellfun(@nanmin,t_full));
t_max=nanmax(cellfun(@nanmax,t_full));
sst_anom=NaN(775,845,length(t_min:t_max));
for i=1:258
    if i<258
        idx_here=(1:3)+(i-1)*3;
    else
        idx_here=772:775;
    end
    mfile=matfile(['/dfs6/pub/zijiz26/nwa/hw3/sst_anom3d' num2str(i)]);
    sst_anom(idx_here,:,:)=mfile.sst_anom(:,:,t_min:t_max);
end

load(['radius_sc' num2str(d)]);
warning off
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
        hw_here=sst_anom(:,:,t_here(j));
        ind_here=sub2ind([775 845],xloc{j},yloc{j});
        img_here=zeros(775,845);
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
            
            if nansum(in | on)>2 && length(unique(lond(in | on)))==1
                sst_b_here=interp1(latd(in | on),hw_here(in | on),y_full(:),'linear','extrap');
                sst_b_here=reshape(sst_b_here,res,res);
            elseif nansum(in | on)>2 && length(unique(latd(in | on)))==1
                sst_b_here=interp1(lond(in | on),hw_here(in | on),x_full(:),'linear','extrap');
                sst_b_here=reshape(sst_b_here,res,res);
            
            elseif nansum(in | on)>2
                F=scatteredInterpolant(lond(in | on),latd(in | on),hw_here(in | on),'linear','none');
                sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
                
            else
                sst_b_here=nanmean(hw_here(in | on));
            end
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
save(['spn_sc' num2str(d)],'sst_sp','-v7.3');

%% radius - stc
% d = depth_s;
% idx_here=(1:5)+(d-1)*5;
% 
% load('areas');
% area_used=areas;
% 
% load('bw3d5r');
% voxellist=regionprops3(bw,'VoxelList');
% objs=voxellist.VoxelList(idx_here,:);
% tracks_new=objs;
% tracks=objs;
% load('lonlat');
% lon_used=xh;
% lat_used=yh;
% 
% radius_sp=cell(length(tracks),1);
% [lat2,lon2]=meshgrid(lat_used,lon_used);
% 
% tmin=nanmin(cellfun(@(x)nanmin(x(:,3)),objs));
% tmax=nanmax(cellfun(@(x)nanmax(x(:,3)),objs));
% radius_s=NaN(length(tracks),1);
% 
% for i=1:length(tracks)
%     %tic
%     %         hw_here=sst_anom(:,:,t_here(j));
%     ind_here=sub2ind([775 845],tracks{i}(:,1),tracks{i}(:,2));
%     img_here=zeros(775,845);
%     img_here(ind_here)=1;
%     bw=bwconncomp(img_here,8);
%     %sst_b=NaN(100,100,bw.NumObjects);
%     radius_b=NaN(1,1);
%     for b=1:1
%         idx_here=ind_here;
%         center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%         %            center_b(b,:)=center_here;
%         dist_here=geodist(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%         
%         %             area_here=length(idx_here)*0.0625;
%         %             area_b(b)=area_here;
%         
%         radius_b(b)=nanmax(dist_here);
%     end
%     %         sst_s(:,:,j)=nansum(sst_b.*reshape(area_b,1,1,length(area_b)),3)./nansum(area_b);
%     %center_s(j,:)=nansum(center_b.*area_b,1)./nansum(area_b);
%     %         x_s{j}=x_b;y_s{j}=y_b;
%     %area_s{j}=area_b;
%     radius_s(i)=radius_b;
%     
%     %     sst_sp{i}=sst_s;
%     %     x_sp{i}=x_s;
%     %     y_sp{i}=y_s;
%     %     area_sp{i}=area_s;
%     %toc
% end
%  save(['radius_stc' num2str(d)],'radius_s','-v7.3');
 
 %% spn stc
d = depth_s;
idx_here=(1:5)+(d-1)*5;
load('areas');
load('lonlat');
area_used=areas;
lon_used=xh;
lat_used=yh;
%addpath /scratch/v45/zz6006/hw_dynamical/hw_new

load('bw3d5r');
voxellist=regionprops3(bw,'VoxelList');
objs=voxellist.VoxelList(idx_here,:);
tracks=objs;

[lat2,lon2]=meshgrid(lat_used,lon_used);
res=50;

t_min=nanmin(cellfun(@(x)nanmin(x(:,3)),objs));
t_max=nanmax(cellfun(@(x)nanmax(x(:,3)),objs));
sst_anom=NaN(775,845,length(t_min:t_max));
for i=1:258
    if i<258
        idx_here=(1:3)+(i-1)*3;
    else
        idx_here=772:775;
    end
    mfile=matfile(['/dfs6/pub/zijiz26/nwa/hw3/sst_anom3d' num2str(i)]);
    sst_anom(idx_here,:,:)=mfile.sst_anom(:,:,t_min:t_max);
end
load(['radius_stc' num2str(d)]);
warning off
sst_sp=cell(5,1);
tic
for i=1:length(objs)
    t_here=unique(tracks{i}(:,3))-t_min+1;
    sst_s=NaN(res,res,length(t_here));
    center_s=NaN(length(t_here),2);
    xloc=tracks{i}(:,2);
    yloc=tracks{i}(:,1);
    for j=1:length(t_here)
        hw_here=sst_anom(:,:,t_here(j));
        ind_here=sub2ind([775 845],xloc,yloc);
        img_here=zeros(775,845);
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
            
            if nansum(in | on)>2 && length(unique(lond(in | on)))==1
                sst_b_here=interp1(latd(in | on),hw_here(in | on),y_full(:),'linear','extrap');
                sst_b_here=reshape(sst_b_here,res,res);
            elseif nansum(in | on)>2 && length(unique(latd(in | on)))==1
                sst_b_here=interp1(lond(in | on),hw_here(in | on),x_full(:),'linear','extrap');
                sst_b_here=reshape(sst_b_here,res,res);
                
            elseif nansum(in | on)>2
                F=scatteredInterpolant(lond(in | on),latd(in | on),hw_here(in | on),'linear','none');
                sst_b_here=reshape(F(x_full(:),y_full(:)),res,res);
                
            else
                sst_b_here=nanmean(hw_here(in | on));
            end
            %sst_b_here(~in & ~on)=0;
            sst_b(:,:,b)=sst_b_here;
            x_b(:,:,b)=x_full;
            y_b(:,:,b)=y_full;
            
        end
        sst_s(:,:,j)=sst_b;
        center_s(j,:)=center_b;
    end
    sst_sp{i}=sst_s;
end
toc
save(['spn_stc' num2str(d)],'sst_sp','-v7.3');

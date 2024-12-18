%% MHW radius
% for d=1:5
%     load(['/g/data/v45/zz6006/mhw_dynamical/tracksr_' num2str(d)]);
%     load('lonlat_access');
%     tracks_new=squeeze(struct2cell(tracks));
%     lon_used=lon_oi;
%     lat_used=lat_oi;
%     
%     radius_sp=cell(length(tracks),1);
%     [lat2,lon2]=meshgrid(lat_used,lon_used);
%     
%     for i=1:length(tracks)
%         tic
%         t_here=tracks_new{1,i};
%         radius_s=cell(length(t_here),1);
%         xloc=tracks_new{2,i};
%         yloc=tracks_new{3,i};
%         for j=1:length(t_here)
%             %         mhw_here=sst_anom(:,:,t_here(j));
%             ind_here=sub2ind([400 251],xloc{j},yloc{j});
%             img_here=zeros(400,251);
%             img_here(ind_here)=1;
%             bw=bwconncomp(img_here,8);
%             %sst_b=NaN(100,100,bw.NumObjects);
%             radius_b=NaN(1,1);
%             for b=1:1
%                 idx_here=ind_here;
%                 center_here=[nanmean(lon2(idx_here)) nanmean(lat2(idx_here))];
%                 %            center_b(b,:)=center_here;
%                 dist_here=geodist(lat2(idx_here),center_here(2),lon2(idx_here),center_here(1));
%                 
%                 %             area_here=length(idx_here)*0.0625;
%                 %             area_b(b)=area_here;
%                 
%                 radius_b(b)=nanmax(dist_here);
%             end
%             %         sst_s(:,:,j)=nansum(sst_b.*reshape(area_b,1,1,length(area_b)),3)./nansum(area_b);
%             %center_s(j,:)=nansum(center_b.*area_b,1)./nansum(area_b);
%             %         x_s{j}=x_b;y_s{j}=y_b;
%             %area_s{j}=area_b;
%             radius_s{j}=radius_b;
%         end
%         %     sst_sp{i}=sst_s;
%         %     x_sp{i}=x_s;
%         %     y_sp{i}=y_s;
%         %     area_sp{i}=area_s;
%         radius_sp{i}=radius_s;
%         toc
%     end
%     radius_s=NaN(length(radius_sp),1);
%     for i=1:length(radius_sp)
%         radius_here=radius_sp{i};
%         radius_here=cellfun(@nanmax,radius_here);
%         radius_s(i)=nanmax(radius_here);
%     end
%     save(['/g/data/v45/zz6006/mhw_dynamical/radiusr_' num2str(d)],'radius_s');
% end

%% spn
d = depth_s;
load(['/g/data/v45/zz6006/mhw_dynamical/radiusr_' num2str(d)]);
p='/g/data/v45/zz6006/mhw_dynamical/';
load(['/g/data/v45/zz6006/mhw_dynamical/tracksr_' num2str(d)]);
tracks_new=squeeze(struct2cell(tracks));
load('lonlat_access');
lon_used=lon_oi;
lat_used=lat_oi;
%load([p 'MHW_access'],'mclim');
%load(['/g/data/v45/zz6006/mhw_dynamical/sst_full_access']);
load('adv_anom');
idx_end=datenum(2018,12,31)-datenum(1987,1,1)+1;
% sst_anom=adv_anom;
date_used=datevec(datenum(1987,1,1):datenum(2018,12,31));
date_used(:,1)=2000;
dayofyear=day(datetime(date_used),'dayofyear');
%sst_anom=sst_full-mclim(:,:,dayofyear);
sst_sp=cell(length(tracks),1);
x_sp=cell(length(tracks),1);
y_sp=cell(length(tracks),1);
area_sp=cell(length(tracks),1);
[lat2,lon2]=meshgrid(lat_used,lon_used);
res=50;

for i=1:length(tracks)
    tic
    t_here=tracks_new{1,i};
    sst_s=NaN(res,res,length(t_here));
    x_s=cell(length(t_here),1);
    y_s=cell(length(t_here),1);
    center_s=NaN(length(t_here),2);
    area_s=cell(length(t_here),1);
    xloc=tracks_new{2,i};
    yloc=tracks_new{3,i};
    for j=1:length(t_here)
        mhw_here=sst_anom(:,:,t_here(j));
        ind_here=sub2ind([400 251],xloc{j},yloc{j});
        img_here=zeros(400,251);
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
    toc
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

 save(['/g/data/v45/zz6006/mhw_dynamical/advspnr_' num2str(d)],'sst_spn','-v7.3');
 
 %% MHW advection classifier
d = depth_s;
load('hb_anom');
dTdt_anom=smoothdata(dTdt_anom,3,'movmean',11);
dTdt_anom=smoothdata(dTdt_anom,3,'movmean',11);
udTdx_anom=smoothdata(udTdx_anom,3,'movmean',11)*3600*24;
vdTdy_anom=smoothdata(vdTdy_anom,3,'movmean',11)*3600*24;
adv_anom=-udTdx_anom-vdTdy_anom;
area_full=repmat(area_used,1,1,11688);

prop_adv_full=cell(5,1);
for d=1:5
    tic
    load(['/g/data/v45/zz6006/mhw_dynamical/tracksr_' num2str(d)]);
    prop_adv=NaN(length(tracks),1);
    
    for m=1:length(tracks)
        xloc=tracks(m).xloc;
        yloc=tracks(m).yloc;
        ts=tracks(m).day;
        ts_full=cell(length(ts),1);
        
        for i=1:length(xloc)
            xloc_here=xloc{i};
            ts_here=ts(i)*ones(length(xloc_here),1);
            ts_full{i}=ts_here;
        end
        xloc=cellfun(@double,xloc,'UniformOutput',false);
        yloc=cellfun(@double,yloc,'UniformOutput',false);
        
        xloc=cell2mat(xloc);
        yloc=cell2mat(yloc);
        ts_full=cell2mat(ts_full);
        
        ind=sub2ind([size(dTdt_anom,1),size(dTdt_anom,2) size(dTdt_anom,3)],double(xloc),double(yloc),ts_full);
        dTdt_here=dTdt_anom(ind);
        adv_here=adv_anom(ind);
        res_here=dTdt_here-adv_here;
        idx_here=(dTdt_here.*adv_here)>=0 & abs(adv_here)>=abs(res_here);
        prop_adv(m)=nansum(area_full(idx_here))./nansum(area_full(ind));
    end
    prop_adv_full{d}=prop_adv;
    toc
end

prop_dom=NaN(5,1);

for i=1:5
    prop_dom(i)=nansum(prop_adv_full{i}>0.5)/length(prop_adv_full{i});
end


 
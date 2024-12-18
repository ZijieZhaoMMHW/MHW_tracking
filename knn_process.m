% d= depth_s;
% idx_here=(1:487)+(d-1)*487;
% mfile=matfile(['/g/data/v45/zz6006/mhw_dynamical/MHW_access']);
% mhw_ts=mfile.mhw_ts(:,:,idx_here);
% load('ocean_index');
% 
% res_used=[0.5 1 1.5 2 2.5];
% for i=1:5
%     tic
%     mhw_re=mhwknn(mhw_ts,res_used(i),ocean_index);
%     save(['/g/data/v45/zz6006/mhw_dynamical/mhw_knn' num2str(d) '_' num2str(i)],'mhw_re','-v7.3');
%     toc
% end

load('lonlat');
d = depth_s;
mhwknn=NaN(400,251,11688);
for i=1:24
    load(['/g/data/v45/zz6006/mhw_dynamical/mhw_knn' num2str(i) '_' num2str(d)]);
    idx_here=(1:487)+(i-1)*487;
    mhwknn(:,:,idx_here)=mhw_re;
end

tracks=mhwtrack(mhwknn,lon_access,lat_access,0.1,d*0.5);
save(['/g/data/v45/zz6006/mhw_dynamical/tracksr_' num2str(d)],'tracks','-v7.3');



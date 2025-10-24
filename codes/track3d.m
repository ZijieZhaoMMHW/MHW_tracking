data_sub=NaN(775,845,3637);

for i=1:258
    if i<258
        loc_here=(1:3)+(i-1)*3;
    else
        loc_here=772:775;
    end
    load(['/dfs6/pub/zijiz26/nwa/mhw3/mhw' num2str(i)]);
    data_sub(loc_here,:,:)=mhw_ts;
end

bw=bwconncomp(~isnan(data_sub),26);
save('bwraw','bw','-v7.3');

% load('lonlat');
% xh(xh<=0)=xh(xh<=0)+360;
% tracks=mhwtrack_nouniform(double(~isnan(data_sub)),xh,yh);
% save(['tracks3draw.mat'],'tracks','-v7.3');

% % mhw_full=NaN(775,845,3637);
% % for d=1:363
% %     if d<363
% %         idx_here=(1:10)+(d-1)*10;
% %     else
% %         idx_here=3621:3637;
% %     end
% %     load(['/dfs6/pub/zijiz26/nwa/mhw3/mhwknn3d3' num2str(d)]);
% %     mhw_full(:,:,idx_here)=data_knn;
% % end
% % 
% % mhw_full=double(mhw_full>=0.5);
% % 
% % for i=1:size(mhw_full,3)
% %     mhw_here=mhw_full(:,:,i);
% %     bw=bwconncomp(mhw_here,8);
% %     for j=1:bw.NumObjects
% %         idx_here=bw.PixelIdxList{j};
% %         if length(idx_here)<961
% %             mhw_here(idx_here)=0;
% %         end
% %     end
% %     mhw_full(:,:,i)=mhw_here;
% % end
% % 
% % mhw_full=NaN(775,845,3637);
% % for d=1:363
% %     if d<363
% %         idx_here=(1:10)+(d-1)*10;
% %     else
% %         idx_here=3621:3637;
% %     end
% %     load(['/dfs6/pub/zijiz26/nwa/mhw3/mhwknn3d3mld' num2str(d)]);
% %     mhw_full(:,:,idx_here)=data_knn;
% % end
% % 
% % bw=bwconncomp(double(mhw_full>=0.5),26);
% 
% % mhw_full=NaN(775,845,3637);
% % for d=1:363
% %     if d<363
% %         idx_here=(1:10)+(d-1)*10;
% %     else
% %         idx_here=3621:3637;
% %     end
% %     load(['/dfs6/pub/zijiz26/nwa/mhw3/mhwknn3d5r' num2str(d)]);
% %     mhw_full(:,:,idx_here)=data_knn;
% % end
% % 
% % mhw_full=double(mhw_full>=0.5);
% % 
% % for i=1:size(mhw_full,3)
% %     mhw_here=mhw_full(:,:,i);
% %     bw=bwconncomp(mhw_here,8);
% %     for j=1:bw.NumObjects
% %         idx_here=bw.PixelIdxList{j};
% %         if length(idx_here)<2601
% %             mhw_here(idx_here)=0;
% %         end
% %     end
% %     mhw_full(:,:,i)=mhw_here;
% % end
% % bw=bwconncomp(double(mhw_full>=0.5),26);
% % save(['bw3d5r.mat'],'bw','-v7.3');
% 
% mhw_full=NaN(775,845,3637);
% for d=1:363
%     if d<363
%         idx_here=(1:10)+(d-1)*10;
%     else
%         idx_here=3621:3637;
%     end
%     load(['/dfs6/pub/zijiz26/nwa/mhw3/mhwknn3d3_' num2str(d)]);
%     mhw_full(:,:,idx_here)=data_knn;
% end
% 
% mhw_full=double(mhw_full>=0.5);
% 
% for i=1:size(mhw_full,3)
%     mhw_here=mhw_full(:,:,i);
%     bw=bwconncomp(mhw_here,8);
%     for j=1:bw.NumObjects
%         idx_here=bw.PixelIdxList{j};
%         if length(idx_here)<2601
%             mhw_here(idx_here)=0;
%         end
%     end
%     mhw_full(:,:,i)=mhw_here;
% end
% numse=NaN(size(mhw_full,1),size(mhw_full,2));
% for i=1:size(mhw_full,1)
%     tic
%     for j=1:size(mhw_full,2)
%         ts_here=squeeze(mhw_full(i,j,:));
%         bw=bwconncomp(double(ts_here==1));
%         numse(i,j)=bw.NumObjects;
%     end
%     toc
% end
% save(['numse3'],'numse','-v7.3');
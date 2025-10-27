function tracks=hwtrack_nouniform(hw,lon_used,lat_used,nums)
% Originally written by Sun Di
% Modified and rewritten by Zijie Zhao

HWs = struct('day',{}, 'xloc',{},'yloc',{}); 

%mask_daily = hw;
for i = 1:length(hw(1,1,:))
    count = 0;
    mask = hw(:,:,i);
    D = bwconncomp(mask,8);
    for j = 1:D.NumObjects
        
        [x,y] = ind2sub([size(hw,1),size(hw,2)],D.PixelIdxList{j});
        if length(x)>=nums
        count = count + 1;
        HWs(i).yloc{count,1} = single(y);
        HWs(i).xloc{count,1} = single(x);
        HWs(i).day = i;
        end
    end
    disp(i)
end

for i = 1:length(HWs)
    judge_hw(i) = isempty(HWs(i).day);
end
empty_hw = find(judge_hw==1);
for i = 1:length(empty_hw)
    HWs(empty_hw(i)).day = empty_hw(i);
end

% Input the parameter
% -------------------------------------------------------------------------
para_alpha = 0.5; 
% -------------------------------------------------------------------------

cut_off = 5;
search = struct('day',{}, 'xloc',{},'yloc',{},'ori_day',{},'ori_order',{},'split_num',{},'split_day',{}); 
tracks = struct('day',{}, 'xloc',{},'yloc',{},'ori_day',{},'ori_order',{},'split_num',{},'split_day',{}); 

pi180 = pi/180;
earth_radius = 6378.137;

% ------------------------------ Beginning --------------------------------
for i = 1:length(HWs)
    day = HWs(i).day;
    hw_xloc = HWs(i).xloc;
    hw_yloc = HWs(i).yloc;
    
    % if the first day, open tracks
    count_order = 0;
    if i == 1
        for i2 = 1:length(hw_xloc)
            count_order = count_order + 1;
            search(i2).day = day;
            search(i2).xloc = hw_xloc(i2);
            search(i2).yloc = hw_yloc(i2);
            
            search(i2).ori_day = day;
            search(i2).ori_order = count_order;
            search(i2).split_num = [];
            search(i2).split_day = [];
        end
    else
        % if not the first day, put the xloc and yloc into loc_now
        loc_now = [];
        for i2 = 1:length(hw_xloc)
            loc_now(i2).xloc = hw_xloc(i2);
            loc_now(i2).yloc = hw_yloc(i2);
            
            loc_now(i2).ori_day = i;
            loc_now(i2).ori_order = i2;
        end
        
        count = zeros(length(hw_xloc),1); % for counting times that idx is used
        
        for i2 = 1:length(search)
            % find xloc and yloc in open tracks from the previous day 
            if (day - search(i2).day(end))==1
                
                loc_old(1) = search(i2).xloc(end);
                loc_old(2) = search(i2).yloc(end);
            
                overlap = zeros(size(loc_now));
                judge1 = zeros(size(hw,1),size(hw,2));
                loc_old_x = cell2mat(loc_old(1));
                loc_old_y = cell2mat(loc_old(2));
                for k = 1:length(loc_old_x)
                    judge1(loc_old_x(k),loc_old_y(k)) = 1;
                end

                % loop for all HWs at t2
                for i3 = 1:length(loc_now)
                    judge2 = zeros(size(hw,1),size(hw,2));
                    loc_now_x = cell2mat(loc_now(i3).xloc);
                    loc_now_y = cell2mat(loc_now(i3).yloc);
                    for k2 = 1:length(loc_now_x)
                        judge2(loc_now_x(k2),loc_now_y(k2)) = 1;
                    end
                    
                    % ------------------- Important Part ! ----------------
                    % compute overlap with a fixed t1 when loop for t2
                    overlap(i3) = length(find( judge1 == 1 & judge2 == 1 )) / min(sum(sum(judge2)),sum(sum(judge1)));
                    % -----------------------------------------------------------------------------------------------
                end

                idx = find(overlap >= para_alpha);
                % =========================================================
                % maybe find more than one HWs at t2 overlap>0.5 
                % for splitting
                if ~isempty(idx)
                    if length(idx) > 1
                        search(i2).day=[search(i2).day; day];
                        search(i2).xloc=[search(i2).xloc; hw_xloc(idx(1))];
                        search(i2).yloc=[search(i2).yloc; hw_yloc(idx(1))]; 
                        disp('split!!')
                        disp(i)
                        for i4 = 2:length(idx)
                            search(i2).xloc{end,1} = [search(i2).xloc{end,1}; cell2mat(hw_xloc(idx(i4)))];
                            search(i2).yloc{end,1} = [search(i2).yloc{end,1}; cell2mat(hw_yloc(idx(i4)))];
                        end
                        search(i2).split_num = [search(i2).split_num; length(idx)];
                        search(i2).split_day = [search(i2).split_day; day];
                    else
                        search(i2).day=[search(i2).day; day];
                        search(i2).xloc=[search(i2).xloc; hw_xloc(idx)];
                        search(i2).yloc=[search(i2).yloc; hw_yloc(idx)];
                    end
                end
                % ===================== End splitting ==================
                
                % =========================================================
                % for merging
                count(idx) = count(idx) + 1; 
            end
        end
        
        % ======================== Begin merging ===========================
        % find pervious tracks
        if ~isempty(find(count>1))
            idx_now = find(count>1);
            for m = 1:length(idx_now) % m is for how many HWs in t2 are overlapped with t1
                % find t2 hw in which t1 tracks before
                
                c = 0; % for counting this hw in how many t1 tracks
                loc_now_xloc = cell2mat(loc_now(idx_now(m)).xloc);
                loc_now_yloc = cell2mat(loc_now(idx_now(m)).yloc);
                for n = 1:length(search)
                    search_xloc = cell2mat(search(n).xloc(end));
                    search_yloc = cell2mat(search(n).yloc(end));
                    
                    l = zeros(size(hw,1),size(hw,2));
                    for j5 = 1:length(loc_now_xloc)
                        l(loc_now_xloc(j5),loc_now_yloc(j5))=1;
                    end
                    ll = zeros(size(hw,1),size(hw,2));
                    for j5 = 1:length(search_xloc)
                        ll(search_xloc(j5),search_yloc(j5))=1;
                    end
                    
                    if day == search(n).day(end) &&...
                            ( ...
                            ( length(loc_now_xloc) == length(search_xloc) && isequal(loc_now_xloc, search_xloc) ) || ...
                            ( sum(sum(find(l==1 & ll==1))) == sum(sum(find(l==1))) )...
                             ) % A & (B | C)
                             % for both splitting and merging
                        c = c + 1;
                        old(m,c) = n;
                        
%                         % -------------- for testing -----------------
%                         if day == search(n).day(end) &&...
%                             ( ...
%                             ( length(loc_now_xloc) ~= length(search_xloc)  ) && ...
%                             ( sum(sum(find(l==1 & ll==1))) == sum(sum(find(l==1))) )...
%                              )
%                             disp('attention!!!')
%                             disp(n)
%                         end
%                         % --------------------------------------------
                        
                    end
                end
            end
            
            if m>=2
                disp('Please take care!')
                disp(i)
            end
            
            % for number conservation
            for i5 = 1:length(idx_now)
                count_new_fusion = 0;
                new_fusion = [];
                new_fusion_loc = [];
                include_x = loc_now(idx_now(i5)).xloc{1,1};
                include_y = loc_now(idx_now(i5)).yloc{1,1};
                for j = 1:length(find(old(i5,:)~=0)) % the last one of old maybe 0
                    % part I for those overlapped area
                    past_2d = zeros(size(hw,1),size(hw,2));
                    for j3 = 1:length(search(old(i5,j)).xloc{end-1,1})
                        past_2d(search(old(i5,j)).xloc{end-1,1}(j3),search(old(i5,j)).yloc{end-1,1}(j3)) = 1;
                    end
                    now_2d = zeros(size(hw,1),size(hw,2));
                    for j3 = 1:length(search(old(i5,j)).xloc{end,1})
                        now_2d(search(old(i5,j)).xloc{end,1}(j3),search(old(i5,j)).yloc{end,1}(j3)) = 1;
                    end
                    [x3,y3] = ind2sub([size(hw,1),size(hw,2)],find(past_2d==1 & now_2d==1));
                    
                    for j3 = 1:length(x3)
                        test = include_x==x3(j3) & include_y==y3(j3);
                        include_x(test)=-j;
                        include_y(test)=-j;
                    end
                    
                    mean_x(j) = mean(search(old(i5,j)).xloc{end-1,1});
                    mean_y(j) = mean(search(old(i5,j)).yloc{end-1,1});

                end
                
                 % part II for others
                include_x_0 = include_x(include_x>0);
                include_y_0 = include_y(include_x>0);
                
                if ~isempty(include_x_0)
                    for i6 = 1:length(mean_x)
    %                     a(:,i6) = sqrt( (include_x(include_x>0) - mean_x(i6)).^2 ...
    %                               +( include_y(include_x>0) - mean_y(i6)).^2 );
                        %lat2 = ((mean_y(i6)-1) * res + nanmin(lat_used)) * pi180;
                        lat2 = lat_used(floor(mean_y(i6)))+(mean_y(i6)-floor(mean_y(i6))).*(lat_used(ceil(mean_y(i6)))-...
                            lat_used(floor(mean_y(i6))));
                        
                        %lon2 = ((mean_x(i6)-1) * res + nanmin(lon_used)) * pi180;
                        lon2 = lon_used(floor(mean_x(i6)))+(mean_x(i6)-floor(mean_x(i6))).*(lon_used(ceil(mean_x(i6)))-...
                            lon_used(floor(mean_x(i6))));

                        for ii = 1:length(include_x_0)
                            %lat1 = ((include_y_0(ii)-1) * res + nanmin(lat_used)) * pi180;
                            lat1 = lat_used(floor(include_y_0(ii)))+(include_y_0(ii)-floor(include_y_0(ii))).*(lat_used(ceil(include_y_0(ii)))-...
                                lat_used(floor(include_y_0(ii))));
                            %lon1 = ((include_x_0(ii)-1) * res + nanmin(lon_used)) * pi180;
                            lon1 = lon_used(floor(include_x_0(ii)))+(include_x_0(ii)-floor(include_x_0(ii))).*(lon_used(ceil(include_x_0(ii)))-...
                                lon_used(floor(include_x_0(ii))));
                            dlat = lat1 - lat2;
                            dlon = lon1 - lon2;
                            alpha = (sin(dlat/2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon/2)).^2;
                            angles = 2 * atan2( abs(sqrt(alpha)), abs(sqrt(1-alpha)) );
                            a(ii,i6) = earth_radius * angles;
                        end
    %                     disp('This block is OK !!!')
                    end

                    [s,t] = min(a,[],2);
                    include_x(include_x>0) = -t;
                    include_y(include_y>0) = -t;
                end
                
                % append
                for j = 1:length(find(old(i5,:)~=0))
                    if length(search(old(i5,j)).xloc{end,1})==length(include_x)
                        search(old(i5,j)).xloc{end,1} = search(old(i5,j)).xloc{end,1}(include_x==-j);
                        search(old(i5,j)).yloc{end,1} = search(old(i5,j)).yloc{end,1}(include_y==-j);
                    
                    else
                        % The next is for both split and merge e.g. day 285
                        extract1 = zeros(size(hw,1),size(hw,2));
                        for j1 = 1:length(search(old(i5,j)).xloc{end,1})
                            extract1(search(old(i5,j)).xloc{end,1}(j1),search(old(i5,j)).yloc{end,1}(j1))=1;
                        end
                        extract2 = zeros(size(hw,1),size(hw,2));
                        for j1 = 1:length(loc_now(idx_now(i5)).xloc{1,1})
                            extract2(loc_now(idx_now(i5)).xloc{1,1}(j1),loc_now(idx_now(i5)).yloc{1,1}(j1))=1;
                        end
                        [x1,y1] = ind2sub([size(hw,1),size(hw,2)],find(extract1==1 & extract2==0));
                        [x2,y2] = ind2sub([size(hw,1),size(hw,2)],find(extract1==1 & extract2==1));
                        
                        search(old(i5,j)).xloc{end,1} = [x1;x2(include_x==-j)];
                        search(old(i5,j)).yloc{end,1} = [y1;y2(include_y==-j)];
                       
                        count_new_fusion = count_new_fusion+1;
                        new_fusion_loc(count_new_fusion) = j;
                        new_fusion(count_new_fusion).xloc{1,1} =  x2(include_x==-j);
                        new_fusion(count_new_fusion).yloc{1,1} =  y2(include_y==-j);
                    end
                end
                
                % ------------------------------------------------------------------------------------------------------------
                % The following is for those unconnected according to distance
                if isempty(new_fusion_loc)
                    % if there is no new split, loop for all xloc and yloc
                    % to find connectivity
                    % ------------------------------------------------------------------------
                    for j = 1:length(find(old(i5,:)~=0))
                        split_fusion = zeros(size(hw,1),size(hw,2));
                        for k3 = 1:length(search(old(i5,j)).xloc{end,1})
                            split_fusion(search(old(i5,j)).xloc{end,1}(k3),search(old(i5,j)).yloc{end,1}(k3)) = 1;
                        end
                        D = bwconncomp(split_fusion,8); 
                    % -------------------------------------------------------------------------
                    
                        % if find unconnected contours
                        if D.NumObjects > 1
                            for k7 = 1:D.NumObjects
                                b(k7) =  length(D.PixelIdxList{k7});
                            end
                            [p,q] = max(b);
                            for k4 = 1:D.NumObjects
                                if k4~=q
                                    [x,y] = ind2sub([size(hw,1),size(hw,2)],D.PixelIdxList{k4});
                                    for k5 = 1:length(find(old(i5,:)~=0))
                                        if k5~=j
                                            append_x = [search(old(i5,k5)).xloc{end,1}; x];
                                            append_y = [search(old(i5,k5)).yloc{end,1}; y];

                                            raw_search = zeros(size(hw,1),size(hw,2));
                                            for k6 = 1:length(search(old(i5,k5)).xloc{end,1})
                                                raw_search(search(old(i5,k5)).xloc{end,1}(k6),search(old(i5,k5)).yloc{end,1}(k6)) = 1;
                                            end
                                            D_raw = bwconncomp(raw_search,8); 

                                            test_fusion = zeros(size(hw,1),size(hw,2));
                                            for k6 = 1:length(append_x)
                                                test_fusion(append_x(k6),append_y(k6)) = 1;
                                            end
                                            D1 = bwconncomp(test_fusion,8); 

                                            if D1.NumObjects <= D_raw.NumObjects
                                                % put the merged to the new one
                                                search(old(i5,k5)).xloc{end,1} = append_x;
                                                search(old(i5,k5)).yloc{end,1} = append_y;

                                                % extract the x,y from the old
                                                extract3 = zeros(size(hw,1),size(hw,2));
                                                for j2 = 1:length(x)
                                                    extract3(x(j2),y(j2))=1;
                                                end
                                                extract4 = zeros(size(hw,1),size(hw,2));
                                                for j2 = 1:length(search(old(i5,j)).xloc{end,1})
                                                    extract4(search(old(i5,j)).xloc{end,1}(j2),search(old(i5,j)).yloc{end,1}(j2))=1;
                                                end
                                                [x4,y4] = ind2sub([size(hw,1),size(hw,2)],find(extract3==0 & extract4==1));
                                                search(old(i5,j)).xloc{end,1} = x4;
                                                search(old(i5,j)).yloc{end,1} = y4;

                                                break; % if find a track that the rest is connected, then append into it and stop loop 
                                                % i.e. append to the first one
                                            end
                                        end
                                    end
                                end
                            end

                        end
                    % ---------------------------------------------------------------------------------------------------------------------
                    clear b
                    end
                else
                    % if there are new splits, loop for those old xloc and yloc
                    % to find connectivity
                    % ------------------------------------------------------------------------
                    for j = 1:length(find(old(i5,:)~=0))
                        if isempty(find(j==new_fusion_loc,1))
                            split_fusion = zeros(size(hw,1),size(hw,2));
                            for k3 = 1:length(search(old(i5,j)).xloc{end,1})
                                split_fusion(search(old(i5,j)).xloc{end,1}(k3),search(old(i5,j)).yloc{end,1}(k3)) = 1;
                            end
                            D = bwconncomp(split_fusion,8); 
                        else
                            split_fusion = zeros(size(hw,1),size(hw,2));
                            for k3 = 1:length(new_fusion(find(j==new_fusion_loc)).xloc{1,1})
                                split_fusion(new_fusion(find(j==new_fusion_loc)).xloc{1,1}(k3),new_fusion(find(j==new_fusion_loc)).yloc{1,1}(k3)) = 1;
                            end
                            D = bwconncomp(split_fusion,8); 
                        end
                    % ---------------------------------------------------------------------------
                        
                        % if find unconnected contours
                        if D.NumObjects > 1
                            for k7 = 1:D.NumObjects
                                b(k7) =  length(D.PixelIdxList{k7});
                            end
                            [p,q] = max(b);
                            for k4 = 1:D.NumObjects
                                if k4~=q
                                    [x,y] = ind2sub([size(hw,1),size(hw,2)],D.PixelIdxList{k4});
                                    for k5 = 1:length(find(old(i5,:)~=0))
                                        if k5~=j
                                            append_x = [search(old(i5,k5)).xloc{end,1}; x];
                                            append_y = [search(old(i5,k5)).yloc{end,1}; y];

                                            raw_search = zeros(size(hw,1),size(hw,2));
                                            for k6 = 1:length(search(old(i5,k5)).xloc{end,1})
                                                raw_search(search(old(i5,k5)).xloc{end,1}(k6),search(old(i5,k5)).yloc{end,1}(k6)) = 1;
                                            end
                                            D_raw = bwconncomp(raw_search,8); 

                                            test_fusion = zeros(size(hw,1),size(hw,2));
                                            for k6 = 1:length(append_x)
                                                test_fusion(append_x(k6),append_y(k6)) = 1;
                                            end
                                            D1 = bwconncomp(test_fusion,8); 

                                            if D1.NumObjects <= D_raw.NumObjects
                                                % put the merged to the new one
                                                search(old(i5,k5)).xloc{end,1} = append_x;
                                                search(old(i5,k5)).yloc{end,1} = append_y;

                                                % extract the x,y from the old
                                                extract3 = zeros(size(hw,1),size(hw,2));
                                                for j2 = 1:length(x)
                                                    extract3(x(j2),y(j2))=1;
                                                end
                                                extract4 = zeros(size(hw,1),size(hw,2));
                                                for j2 = 1:length(search(old(i5,j)).xloc{end,1})
                                                    extract4(search(old(i5,j)).xloc{end,1}(j2),search(old(i5,j)).yloc{end,1}(j2))=1;
                                                end
                                                [x4,y4] = ind2sub([size(hw,1),size(hw,2)],find(extract3==0 & extract4==1));
                                                search(old(i5,j)).xloc{end,1} = x4;
                                                search(old(i5,j)).yloc{end,1} = y4;

                                                break; % if find a track that the rest is connected, then append into it and stop loop 
                                                % i.e. append to the first one
                                            end
                                        end
                                    end
                                end
                            end

                        end
                        clear b
                    % ---------------------------------------------------------------------------------------------------------------------
                    end
                end
                    
                clear a mean_x mean_y
            end
            
            
            
%             % append
%             for i5 = 1:length(idx_now)
%                 for j = 2:length(find(old(i5,:)~=0)) % the last one of old maybe 0
%                     first_tracks_day = search(old(i5,1)).day;
%                     
%                     for l2 = 1:length(search(old(i5,j)).day)-1 % the last hw is the same
%                         if ~isempty(find(search(old(i5,j)).day(l2) == first_tracks_day))
%                             loc = find(search(old(i5,j)).day(l2) == first_tracks_day);
%                         end
% 
%                         search(old(i5,1)).xloc{loc,1} = [ search(old(i5,1)).xloc{loc,1}; search(old(i5,j)).xloc{l2,1}];
%                         search(old(i5,1)).yloc{loc,1} = [ search(old(i5,1)).yloc{loc,1}; search(old(i5,j)).yloc{l2,1}];
%                     end
%                 end
%             end
%             find_remove_loc = old(1:length(idx_now),2:end);
%             remove_loc = find_remove_loc(find_remove_loc~=0);% prevent the situation where 0 is existed in old
%             search(remove_loc)=[]; % prevent the order of search changing

        end
        % ================================ End Merging ==============================

        % =================================================================
        % for new tracks
        hw_xloc = hw_xloc(find(count==0));
        hw_yloc = hw_yloc(find(count==0));
        if ~isempty(hw_xloc)
            for i7 = 1:length(hw_xloc)
                search(length(search)+1).day = day;
                search(length(search)).xloc = hw_xloc(i7);
                search(length(search)).yloc = hw_yloc(i7);
                
                search(length(search)).ori_day = day;
                for i8 = 1:length(loc_now)
                    iden_loc_now_xloc = cell2mat(loc_now(i8).xloc);
                    if length(cell2mat(hw_xloc(i7))) == length(iden_loc_now_xloc)
                    search(length(search)).ori_order = i8;
                    end
                end
                
            end
        end     
        % =================================================================
        
        % =================================================================
        % for moved tracks
        moved = [];
        for i2=1:length(search)
            % find dissipated tracks
            moved(i2)= search(i2).day(end) <= day-1;
            if moved(i2)
                % move tracks to closed tracks array
                tracks(length(tracks)+1)=search(i2);
            end
        end
        % remove tracks from open track array
        search(moved==1)=[];
        % =================================================================
        
        clear loc_now old
    end
    disp(i)
end

 % Add tracks that are still in the search array at the end of the
% time-series
for i=1:length(search)
    tracks(length(tracks)+1)=search(i);
end

% %%%%%%%%%%%%%%%%%%%% remove tracks shorter than cut_off days %%%%%%%%%%%%%%
% for i=1:length(tracks)
%     short(i)=length(tracks(i).day)<cut_off;
% end
% tracks(short==1)=[];
% short=sum(short);     

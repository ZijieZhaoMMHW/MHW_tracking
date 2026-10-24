function [tracks,daily_assignment_audit] = mhwtrack(mhw,lon_used,lat_used,nums)
%MHWTRACK Track binary marine heatwaves with periodic longitude.
% tracks = mhwtrack(mhw,lon_used,lat_used,nums)
% [tracks,audit] = mhwtrack(...) additionally returns daily ownership checks.
%
% mhw: full real binary array [longitude,latitude,time], already smoothed.
% lon_used,lat_used: finite coordinate vectors matching the first two axes.
% nums: positive integer minimum size of each initial daily component.
% Longitude is periodic; use this implementation for a full zonal domain.
% NaN is rejected, not silently treated as an event. Mask missing cells first.
%
% Output fields match the original interface: day, xloc, yloc, ori_day,
% ori_order, split_num, split_day. Coordinates are 1-based single indices.
% Tracks are returned in closure order, then active tracks at the right edge.
% Tracks lasting fewer than five days are omitted. No smoothing is done here.
%
% Audit columns: expected cells, assigned cells, within-track duplicates,
% between-track duplicates, missing cells, outside cells, active tracks.
% Daily conservation assertions are enabled even when audit is not requested.
%
% Originally written by Sun Di; modified and rewritten by Zijie Zhao.
% Modified standalone release, 2026-10-18, GPL-3.0; see NOTICE and LICENSE.txt.

narginchk(4,4);
assert((isnumeric(mhw) || islogical(mhw)) && isreal(mhw) && ...
    ~issparse(mhw) && ndims(mhw)<=3,'mhw:input', ...
    'mhw must be a full real numeric or logical [lon,lat,time] array.');
validateattributes(lon_used,{'numeric'},{'vector','real','finite','nonempty'},mfilename,'lon_used',2);
validateattributes(lat_used,{'numeric'},{'vector','real','finite','nonempty'},mfilename,'lat_used',3);
validateattributes(nums,{'numeric'},{'scalar','real','finite','integer','positive'},mfilename,'nums',4);
assert(numel(lon_used)==size(mhw,1) && numel(lat_used)==size(mhw,2), ...
    'mhw:coordinates','Coordinate lengths must match the first two dimensions.');

tracks = struct('day',{},'xloc',{},'yloc',{},'ori_day',{}, ...
    'ori_order',{},'split_num',{},'split_day',{});
search = struct('day',{},'xloc',{},'yloc',{},'ori_day',{}, ...
    'ori_order',{},'split_num',{},'split_day',{},'track_id',{});
history = cell(0,1);
next_id = 1;
nt = size(mhw,3);
daily_assignment_audit = [];
if nargout>1, daily_assignment_audit=zeros(nt,7); end

for day=1:nt
    mask=mhw(:,:,day);
    assert(all(mask(:)==0 | mask(:)==1),'mhw:binary', ...
        'Day %d contains values other than 0 or 1 (including NaN/Inf).',day);
    [search,closed,audit]=mhwtrack_step(logical(mask),lon_used,lat_used,nums,search,day);
    if nargout>1, daily_assignment_audit(day,:)=audit; end

    % Closed tracks have already stored their final snapshot on the last day.
    for k=1:numel(closed)
        id=closed(k).track_id;
        if history{id}.count>=5
            tracks(end+1)=restore_history(closed(k),history{id}); %#ok<AGROW>
        end
        history{id}=[];
    end

    for k=1:numel(search)
        if isempty(search(k).track_id)
            id=next_id; next_id=next_id+1;
            search(k).track_id=id;
            history{id,1}=struct('count',0,'day',zeros(16,1), ...
                'xloc',{cell(16,1)},'yloc',{cell(16,1)});
        else
            id=search(k).track_id;
        end
        n=history{id}.count+1;
        if n>numel(history{id}.day)
            capacity=2*numel(history{id}.day);
            history{id}.day(capacity,1)=0;
            history{id}.xloc{capacity,1}=[];
            history{id}.yloc{capacity,1}=[];
        end
        history{id}.day(n)=day;
        history{id}.xloc{n}=search(k).xloc{end};
        history{id}.yloc{n}=search(k).yloc{end};
        history{id}.count=n;

        % Earlier coordinates are output history, not tracking state.
        first=max(1,numel(search(k).day)-1);
        search(k).day=search(k).day(first:end);
        search(k).xloc=search(k).xloc(first:end);
        search(k).yloc=search(k).yloc(first:end);
    end
end

for k=1:numel(search)
    id=search(k).track_id;
    if history{id}.count>=5
        tracks(end+1)=restore_history(search(k),history{id}); %#ok<AGROW>
    end
    history{id}=[];
end
end

function track=restore_history(active,stored)
track=rmfield(active,'track_id');
n=stored.count;
track.day=stored.day(1:n);
track.xloc=stored.xloc(1:n);
track.yloc=stored.yloc(1:n);
end

function [search,closed,audit]=mhwtrack_step(mask,lon_used,lat_used,nums,search,day)
% Generated day logic: same order, thresholds, merge assignment and connectivity.
% History before yesterday is not read by this algorithm and is stored on disk.
mhw=mask; i=day; para_alpha=.5; pi180=pi/180; earth_radius=6378.137;
objects=struct('day',day,'xloc',{{}},'yloc',{{}});
D=bwconncomp_periodic_lon(mask); count=0;
for j=1:D.NumObjects
    [x,y]=ind2sub(size(mask),D.PixelIdxList{j});
    if numel(x)>=nums
        count=count+1; objects.xloc{count,1}=single(x); objects.yloc{count,1}=single(y);
    end
end
mhw_xloc=objects.xloc; mhw_yloc=objects.yloc;
closed=search([]);
    % if the first day, open tracks
    count_order = 0;
    if i == 1
        for i2 = 1:length(mhw_xloc)
            count_order = count_order + 1;
            search(i2).day = day;
            search(i2).xloc = mhw_xloc(i2);
            search(i2).yloc = mhw_yloc(i2);
            
            search(i2).ori_day = day;
            search(i2).ori_order = count_order;
            search(i2).split_num = [];
            search(i2).split_day = [];
        end
    else
        % if not the first day, put the xloc and yloc into loc_now
        loc_now = [];
        for i2 = 1:length(mhw_xloc)
            loc_now(i2).xloc = mhw_xloc(i2);
            loc_now(i2).yloc = mhw_yloc(i2);
            
            loc_now(i2).ori_day = i;
            loc_now(i2).ori_order = i2;
        end
        
        count = zeros(length(mhw_xloc),1); % for counting times that idx is used
        now_labels=zeros(size(mhw,1),size(mhw,2),'uint32');
        now_sizes=zeros(length(loc_now),1);
        for object_index=1:length(loc_now)
            object_cells=sub2ind(size(now_labels),double(cell2mat(loc_now(object_index).xloc)), ...
                double(cell2mat(loc_now(object_index).yloc)));
            now_labels(object_cells)=object_index;
            now_sizes(object_index)=numel(unique(object_cells));
        end
        
        for i2 = 1:length(search)
            % find xloc and yloc in open tracks from the previous day 
            if (day - search(i2).day(end))==1
                
                loc_old(1) = search(i2).xloc(end);
                loc_old(2) = search(i2).yloc(end);
            
                old_cells=unique(sub2ind(size(now_labels),double(cell2mat(loc_old(1))), ...
                    double(cell2mat(loc_old(2)))));
                hits=accumarray(double(now_labels(old_cells))+1,1,[length(loc_now)+1,1]);
                overlap=reshape(hits(2:end)./min(numel(old_cells),now_sizes),size(loc_now));

                idx = find(overlap >= para_alpha);
                % =========================================================
                % maybe find more than one mhws at t2 overlap>0.5 
                % for splitting
                if ~isempty(idx)
                    if length(idx) > 1
                        search(i2).day=[search(i2).day; day];
                        search(i2).xloc=[search(i2).xloc; mhw_xloc(idx(1))];
                        search(i2).yloc=[search(i2).yloc; mhw_yloc(idx(1))]; 
                        for i4 = 2:length(idx)
                            search(i2).xloc{end,1} = [search(i2).xloc{end,1}; cell2mat(mhw_xloc(idx(i4)))];
                            search(i2).yloc{end,1} = [search(i2).yloc{end,1}; cell2mat(mhw_yloc(idx(i4)))];
                        end
                        search(i2).split_num = [search(i2).split_num; length(idx)];
                        search(i2).split_day = [search(i2).split_day; day];
                    else
                        search(i2).day=[search(i2).day; day];
                        search(i2).xloc=[search(i2).xloc; mhw_xloc(idx)];
                        search(i2).yloc=[search(i2).yloc; mhw_yloc(idx)];
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
            for m = 1:length(idx_now) % m is for how many MHWs in t2 are overlapped with t1
                % find t2 MHW in which t1 tracks before
                
                c = 0; % for counting this MHW in how many t1 tracks
                loc_now_xloc = cell2mat(loc_now(idx_now(m)).xloc);
                loc_now_yloc = cell2mat(loc_now(idx_now(m)).yloc);
                for n = 1:length(search)
                    search_xloc = cell2mat(search(n).xloc(end));
                    search_yloc = cell2mat(search(n).yloc(end));
                    
                    l = zeros(size(mhw,1),size(mhw,2));
                    for j5 = 1:length(loc_now_xloc)
                        l(loc_now_xloc(j5),loc_now_yloc(j5))=1;
                    end
                    ll = zeros(size(mhw,1),size(mhw,2));
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
                    past_2d = zeros(size(mhw,1),size(mhw,2));
                    for j3 = 1:length(search(old(i5,j)).xloc{end-1,1})
                        past_2d(search(old(i5,j)).xloc{end-1,1}(j3),search(old(i5,j)).yloc{end-1,1}(j3)) = 1;
                    end
                    now_2d = zeros(size(mhw,1),size(mhw,2));
                    for j3 = 1:length(search(old(i5,j)).xloc{end,1})
                        now_2d(search(old(i5,j)).xloc{end,1}(j3),search(old(i5,j)).yloc{end,1}(j3)) = 1;
                    end
                    [x3,y3] = ind2sub([size(mhw,1),size(mhw,2)],find(past_2d==1 & now_2d==1));
                    
                    available=find(include_x>0);
                    candidate_ids=double(include_x(available))+(double(include_y(available))-1)*size(mhw,1);
                    overlap_ids=double(x3)+(double(y3)-1)*size(mhw,1);
                    claim=available(ismember(candidate_ids,overlap_ids));
                    include_x(claim)=-j;
                    include_y(claim)=-j;
                    
                    mean_x(j) = mean_lon_index_periodic(search(old(i5,j)).xloc{end-1,1},size(mhw,1));
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
                        lat2 = lat2*pi180;
                        
                        %lon2 = ((mean_x(i6)-1) * res + nanmin(lon_used)) * pi180;
                        lon2 = interp_lon_periodic(lon_used,mean_x(i6));
                        lon2=lon2*pi180;

                        for ii = 1:length(include_x_0)
                            %lat1 = ((include_y_0(ii)-1) * res + nanmin(lat_used)) * pi180;
                            lat1 = lat_used(floor(include_y_0(ii)))+(include_y_0(ii)-floor(include_y_0(ii))).*(lat_used(ceil(include_y_0(ii)))-...
                                lat_used(floor(include_y_0(ii))));
                            lat1=lat1*pi180;
                            %lon1 = ((include_x_0(ii)-1) * res + nanmin(lon_used)) * pi180;
                            lon1 = interp_lon_periodic(lon_used,include_x_0(ii));
                            lon1=lon1*pi180;
                            dlat = lat1 - lat2;
                            dlon = wrap_to_pi(lon1 - lon2);
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
                        keep=mhw_owner_mask(search(old(i5,j)).xloc{end,1},search(old(i5,j)).yloc{end,1}, ...
                            loc_now(idx_now(i5)).xloc{1,1},loc_now(idx_now(i5)).yloc{1,1},include_x,j,size(mhw,1));
                        search(old(i5,j)).xloc{end,1} = search(old(i5,j)).xloc{end,1}(keep);
                        search(old(i5,j)).yloc{end,1} = search(old(i5,j)).yloc{end,1}(keep);
                    
                    else
                        % The next is for both split and merge e.g. day 285
                        extract1 = zeros(size(mhw,1),size(mhw,2));
                        for j1 = 1:length(search(old(i5,j)).xloc{end,1})
                            extract1(search(old(i5,j)).xloc{end,1}(j1),search(old(i5,j)).yloc{end,1}(j1))=1;
                        end
                        extract2 = zeros(size(mhw,1),size(mhw,2));
                        for j1 = 1:length(loc_now(idx_now(i5)).xloc{1,1})
                            extract2(loc_now(idx_now(i5)).xloc{1,1}(j1),loc_now(idx_now(i5)).yloc{1,1}(j1))=1;
                        end
                        [x1,y1] = ind2sub([size(mhw,1),size(mhw,2)],find(extract1==1 & extract2==0));
                        [x2,y2] = ind2sub([size(mhw,1),size(mhw,2)],find(extract1==1 & extract2==1));
                        
                        keep=mhw_owner_mask(x2,y2,loc_now(idx_now(i5)).xloc{1,1}, ...
                            loc_now(idx_now(i5)).yloc{1,1},include_x,j,size(mhw,1));
                        search(old(i5,j)).xloc{end,1} = [x1;x2(keep)];
                        search(old(i5,j)).yloc{end,1} = [y1;y2(keep)];
                       
                        count_new_fusion = count_new_fusion+1;
                        new_fusion_loc(count_new_fusion) = j;
                        new_fusion(count_new_fusion).xloc{1,1} = x2(keep);
                        new_fusion(count_new_fusion).yloc{1,1} = y2(keep);
                    end
                end
                
                % ------------------------------------------------------------------------------------------------------------
                % The following is for those unconnected according to distance
                if isempty(new_fusion_loc)
                    % if there is no new split, loop for all xloc and yloc
                    % to find connectivity
                    % ------------------------------------------------------------------------
                    for j = 1:length(find(old(i5,:)~=0))
                        split_fusion = zeros(size(mhw,1),size(mhw,2));
                        for k3 = 1:length(search(old(i5,j)).xloc{end,1})
                            split_fusion(search(old(i5,j)).xloc{end,1}(k3),search(old(i5,j)).yloc{end,1}(k3)) = 1;
                        end
                        D = bwconncomp_periodic_lon(split_fusion); 
                    % -------------------------------------------------------------------------
                    
                        % if find unconnected contours
                        if D.NumObjects > 1
                            for k7 = 1:D.NumObjects
                                b(k7) =  length(D.PixelIdxList{k7});
                            end
                            [p,q] = max(b);
                            for k4 = 1:D.NumObjects
                                if k4~=q
                                    [x,y] = ind2sub([size(mhw,1),size(mhw,2)],D.PixelIdxList{k4});
                                    for k5 = 1:length(find(old(i5,:)~=0))
                                        if k5~=j
                                            append_x = [search(old(i5,k5)).xloc{end,1}; x];
                                            append_y = [search(old(i5,k5)).yloc{end,1}; y];

                                            raw_search = zeros(size(mhw,1),size(mhw,2));
                                            for k6 = 1:length(search(old(i5,k5)).xloc{end,1})
                                                raw_search(search(old(i5,k5)).xloc{end,1}(k6),search(old(i5,k5)).yloc{end,1}(k6)) = 1;
                                            end
                                            D_raw = bwconncomp_periodic_lon(raw_search); 

                                            test_fusion = zeros(size(mhw,1),size(mhw,2));
                                            for k6 = 1:length(append_x)
                                                test_fusion(append_x(k6),append_y(k6)) = 1;
                                            end
                                            D1 = bwconncomp_periodic_lon(test_fusion); 

                                            if D1.NumObjects <= D_raw.NumObjects
                                                % put the merged to the new one
                                                search(old(i5,k5)).xloc{end,1} = append_x;
                                                search(old(i5,k5)).yloc{end,1} = append_y;

                                                % extract the x,y from the old
                                                extract3 = zeros(size(mhw,1),size(mhw,2));
                                                for j2 = 1:length(x)
                                                    extract3(x(j2),y(j2))=1;
                                                end
                                                extract4 = zeros(size(mhw,1),size(mhw,2));
                                                for j2 = 1:length(search(old(i5,j)).xloc{end,1})
                                                    extract4(search(old(i5,j)).xloc{end,1}(j2),search(old(i5,j)).yloc{end,1}(j2))=1;
                                                end
                                                [x4,y4] = ind2sub([size(mhw,1),size(mhw,2)],find(extract3==0 & extract4==1));
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
                            split_fusion = zeros(size(mhw,1),size(mhw,2));
                            for k3 = 1:length(search(old(i5,j)).xloc{end,1})
                                split_fusion(search(old(i5,j)).xloc{end,1}(k3),search(old(i5,j)).yloc{end,1}(k3)) = 1;
                            end
                            D = bwconncomp_periodic_lon(split_fusion); 
                        else
                            split_fusion = zeros(size(mhw,1),size(mhw,2));
                            for k3 = 1:length(new_fusion(find(j==new_fusion_loc)).xloc{1,1})
                                split_fusion(new_fusion(find(j==new_fusion_loc)).xloc{1,1}(k3),new_fusion(find(j==new_fusion_loc)).yloc{1,1}(k3)) = 1;
                            end
                            D = bwconncomp_periodic_lon(split_fusion); 
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
                                    [x,y] = ind2sub([size(mhw,1),size(mhw,2)],D.PixelIdxList{k4});
                                    for k5 = 1:length(find(old(i5,:)~=0))
                                        if k5~=j
                                            append_x = [search(old(i5,k5)).xloc{end,1}; x];
                                            append_y = [search(old(i5,k5)).yloc{end,1}; y];

                                            raw_search = zeros(size(mhw,1),size(mhw,2));
                                            for k6 = 1:length(search(old(i5,k5)).xloc{end,1})
                                                raw_search(search(old(i5,k5)).xloc{end,1}(k6),search(old(i5,k5)).yloc{end,1}(k6)) = 1;
                                            end
                                            D_raw = bwconncomp_periodic_lon(raw_search); 

                                            test_fusion = zeros(size(mhw,1),size(mhw,2));
                                            for k6 = 1:length(append_x)
                                                test_fusion(append_x(k6),append_y(k6)) = 1;
                                            end
                                            D1 = bwconncomp_periodic_lon(test_fusion); 

                                            if D1.NumObjects <= D_raw.NumObjects
                                                % put the merged to the new one
                                                search(old(i5,k5)).xloc{end,1} = append_x;
                                                search(old(i5,k5)).yloc{end,1} = append_y;

                                                % extract the x,y from the old
                                                extract3 = zeros(size(mhw,1),size(mhw,2));
                                                for j2 = 1:length(x)
                                                    extract3(x(j2),y(j2))=1;
                                                end
                                                extract4 = zeros(size(mhw,1),size(mhw,2));
                                                for j2 = 1:length(search(old(i5,j)).xloc{end,1})
                                                    extract4(search(old(i5,j)).xloc{end,1}(j2),search(old(i5,j)).yloc{end,1}(j2))=1;
                                                end
                                                [x4,y4] = ind2sub([size(mhw,1),size(mhw,2)],find(extract3==0 & extract4==1));
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
%                     for l2 = 1:length(search(old(i5,j)).day)-1 % the last MHW is the same
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
        mhw_xloc = mhw_xloc(find(count==0));
        mhw_yloc = mhw_yloc(find(count==0));
        if ~isempty(mhw_xloc)
            for i7 = 1:length(mhw_xloc)
                search(length(search)+1).day = day;
                search(length(search)).xloc = mhw_xloc(i7);
                search(length(search)).yloc = mhw_yloc(i7);
                
                search(length(search)).ori_day = day;
                for i8 = 1:length(loc_now)
                    iden_loc_now_xloc = cell2mat(loc_now(i8).xloc);
                    if length(cell2mat(mhw_xloc(i7))) == length(iden_loc_now_xloc)
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
                closed(length(closed)+1)=search(i2);
            end
        end
        % remove tracks from open track array
        search(moved==1)=[];
        % =================================================================
        
        clear loc_now old
    end
audit=mhw_audit_day(search,objects,day,size(mask,1),size(mask,2));
end

function D = bwconncomp_periodic_lon(mask)
% 8-connected components with periodic wrapping along the longitude axis
% (the first dimension of mhw: lon x lat x time).
D = bwconncomp(mask,8);
nx = size(mask,1);
ny = size(mask,2);

if nx < 2 || D.NumObjects < 2
    return
end

labels = zeros(size(mask),'uint32');
for obj = 1:D.NumObjects
    labels(D.PixelIdxList{obj}) = obj;
end

parent = 1:D.NumObjects;
did_merge = false;
for y = 1:ny
    label1 = double(labels(1,y));
    if label1 == 0
        continue
    end

    y_range = max(1,y-1):min(ny,y+1);
    edge_labels = double(labels(nx,y_range));
    edge_labels = edge_labels(edge_labels > 0);
    for k = 1:numel(edge_labels)
        label2 = edge_labels(k);
        root1 = find_root(parent,label1);
        root2 = find_root(parent,label2);
        if root1 ~= root2
            parent(root2) = root1;
            did_merge = true;
        end
    end
end

if ~did_merge
    return
end

root_map = zeros(1,D.NumObjects);
new_pixel_idx_list = {};
for obj = 1:D.NumObjects
    root = find_root(parent,obj);
    if root_map(root) == 0
        new_pixel_idx_list{end+1,1} = D.PixelIdxList{obj}; %#ok<AGROW>
        root_map(root) = numel(new_pixel_idx_list);
    else
        new_idx = root_map(root);
        new_pixel_idx_list{new_idx,1} = [new_pixel_idx_list{new_idx,1}; D.PixelIdxList{obj}]; %#ok<AGROW>
    end
end

D.PixelIdxList = new_pixel_idx_list;
D.NumObjects = numel(new_pixel_idx_list);
end

function root = find_root(parent,idx)
root = idx;
while parent(root) ~= root
    root = parent(root);
end
end

function idx_mean = mean_lon_index_periodic(idx,nx)
idx = double(idx(:));
theta = 2*pi*(idx-1)./nx;
x = mean(cos(theta));
y = mean(sin(theta));

if x == 0 && y == 0
    idx_mean = mean(idx);
    return
end

theta_mean = atan2(y,x);
if theta_mean < 0
    theta_mean = theta_mean + 2*pi;
end

idx_mean = 1 + theta_mean/(2*pi)*nx;
if idx_mean > nx
    idx_mean = idx_mean - nx;
end
end

function lon = interp_lon_periodic(lon_used,idx)
lon_used = double(lon_used(:));
nx = numel(lon_used);
idx_floor = floor(double(idx));
frac = double(idx) - idx_floor;

i1 = mod(idx_floor-1,nx) + 1;
i2 = mod(idx_floor,nx) + 1;
lon1 = lon_used(i1);
lon2 = lon_used(i2);

dlon = lon2 - lon1;
if dlon > 180
    lon2 = lon2 - 360;
elseif dlon < -180
    lon2 = lon2 + 360;
end

lon = lon1 + frac*(lon2-lon1);
end

function angle = wrap_to_pi(angle)
angle = mod(angle + pi,2*pi) - pi;
end

function keep=mhw_owner_mask(query_x,query_y,component_x,component_y,owner,j,nx)
% Ownership follows a grid-cell ID, never the incidental coordinate ordering.
assigned=double(component_x(owner==-j))+(double(component_y(owner==-j))-1)*nx;
query=double(query_x)+(double(query_y)-1)*nx;
keep=ismember(query,assigned);
end

function record=mhw_audit_day(search,objects,day,nx,ny)
% Test conservation before the final lifetime filter; do not repair the data.
expected=false(nx,ny); seen=zeros(nx,ny,'uint16');
for k=1:numel(objects.xloc)
    ids=double(objects.xloc{k})+(double(objects.yloc{k})-1)*nx;
    expected(ids)=true;
end
within=0; between=0; active=0; total=0;
for k=1:numel(search)
    if search(k).day(end)~=day, continue; end
    active=active+1;
    ids=double(search(k).xloc{end})+(double(search(k).yloc{end})-1)*nx;
    unique_ids=unique(ids); within=within+numel(ids)-numel(unique_ids);
    between=between+nnz(seen(unique_ids)); seen(unique_ids)=k;
    total=total+numel(ids);
end
missing=nnz(expected & seen==0); outside=nnz(~expected & seen>0);
record=[nnz(expected) total within between missing outside active];
assert(within==0 && between==0 && missing==0 && outside==0, ...
    'mhw:ownership','Day %d: within=%d across=%d missing=%d outside=%d', ...
    day,within,between,missing,outside);
end

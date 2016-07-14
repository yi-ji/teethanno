function [F,X,frag_color] = annotate_teeth(filepath,anchor_points)

%anchor_points = [-24.71 -44.14 -4.80;-5.30 -2.78 -3.57;30.3 -43.01 -7.985];%lower_left
%anchor_points = [30.7 -37.64 -0.81;-33.35 -36.14 -1.93;0.243 6.627 -5.015];%lower_right
%anchor_points = [-28.22 -36.63 -12.85;0.302 2.174 -4.3;29.5 -35.72 -14.8];%upper_right

%filepath = 'lower_right.obj';

[F,X,status] = readSmf(filepath);
out=tricurv_v01(F,X);

plane = [anchor_points(:,1:2) [1 1 1]'] \ anchor_points(:,3);
plane_dist = abs(X(:,1)*plane(1)+X(:,2)*plane(2)+plane(3)-X(:,3));
plane_dist_threshold = max(plane_dist)*0.3;

flag = plane_dist<plane_dist_threshold;
X = X(flag,:);
out.km = out.km(flag,:);
F = F(flag(F(:,1))==1 & flag(F(:,2))==1 & flag(F(:,3))==1, :);
F_idx = zeros(1,length(flag));
cnt = 0;
for i = 1:length(flag)
    if flag(i) == 1
        cnt = cnt + 1;
        F_idx(i) = cnt;
    end
end

for i = 1:length(F)
    F(i,1) = F_idx(F(i,1));
    F(i,2) = F_idx(F(i,2));
    F(i,3) = F_idx(F(i,3));
end

fprintf('Object cut done.\n');

dist_average = 0;
for i = 1:4:length(F)
    dist_average = dist_average + norm(X(F(i,1),:)-X(F(i,2),:));
    dist_average = dist_average + norm(X(F(i,2),:)-X(F(i,3),:));
    dist_average = dist_average + norm(X(F(i,1),:)-X(F(i,3),:));
end
dist_average = dist_average / (length(F)*3/4);

curv_threshold = -0.1;
curv_upper_bound = -0.05;
dist_threshold = 2.2*dist_average;
step_threshold = length(X)*0.1;
demarcation_ratio = 80;

nodes = (1:length(out.km)) .* (out.km < curv_threshold)';
nodes(nodes==0) = [];
nodes_num = length(nodes);
center = sum(X(nodes,:))/nodes_num;

[graph,graph_size] = buildE(nodes_num,X(nodes,:)',dist_threshold);
[V,V_num] = selectV(graph,graph_size,nodes,nodes_num,step_threshold);

s=0;t=0;
max_dist=0;
for i = 1:length(V)
    dist=norm(X(V(i),:)-center);
    if dist > max_dist
        s=V(i);
        max_dist=dist;
    end
end

max_dist=0;
for i = 1:length(V)
    if (X(V(i),1)-center(1))*(X(s,1)-center(1)) < 0
        dist=norm(X(V(i),:)-center);
        if dist > max_dist
            t=V(i);
            max_dist=dist;
        end
    end
end

center=(center*2+X(s,:)+X(t,:))/4;
temp_center = center + 20*(center - (X(s,:)+X(t,:))/2);

s_=0;t_=0;
max_dist=0;
for i = 1:length(V)
    dist=norm(X(V(i),:)-temp_center);
    if dist > max_dist
        s_=V(i);
        max_dist=dist;
    end
end

max_dist=0;
for i = 1:length(V)
    if (X(V(i),1)-center(1))*(X(s_,1)-temp_center(1)) < 0
        dist=norm(X(V(i),:)-temp_center);
        if dist > max_dist
            t_=V(i);
            max_dist=dist;
        end
    end
end

std_vec = center - (X(s,:)+X(t,:))/2;
norm_std_vec = norm(std_vec);

fprintf('V-set and point s&t selection done.\n');

cos_range = 1 - sum(std_vec.*(X(s,:)-center))/(norm(X(s,:)-center)*norm_std_vec);
seg_num = 30;
seg_dist = zeros(1,seg_num+1);
seg_points_count = zeros(1,seg_num+1);

for i = 1:V_num
    vec = X(V(i),:) - center;
    cos_val = sum(std_vec.*vec)/(norm(vec)*norm_std_vec);
    temp_cos = 1.0;
    for j = 1:seg_num+1
        if abs(temp_cos - cos_val) < 0.5*cos_range/seg_num
            seg_dist(j) = seg_dist(j) + norm(vec);
            seg_points_count(j) = seg_points_count(j) + 1;
            break;
        end
        temp_cos = temp_cos - cos_range/seg_num;
    end
end

for i = 1:seg_num+1
    seg_dist(i) = seg_dist(i) / seg_points_count(i);
end

cos2dist = csaps(1:-cos_range/seg_num:1-cos_range,seg_dist);

fprintf('Border generation done.\n');

center_dist = zeros(1,length(X));
center_cos = zeros(1,length(X));
for i = 1:length(X)
    vec = X(i,:)-center;
    center_dist(i) = norm(vec);
    center_cos(i) = sum(vec.*std_vec)/(center_dist(i)*norm_std_vec);
end

out.km(out.km>=curv_upper_bound) = curv_upper_bound;
F = int32(F);
border = fnval(cos2dist,center_cos);

[graph,graph_size] = build_graph(F',out.km,center_dist,border,1);
[d,pre] = dijsktra(length(X),graph,graph_size,s_);
inner_nodes = dijsktra([],pre,t_);

[graph,graph_size] = build_graph(F',out.km,center_dist,border,2);
[d,pre] = dijsktra(length(X),graph,graph_size,s_);
outer_nodes = dijsktra([],pre,t_);

fprintf('Contour generation done.\n');

angle_bound = acos(sum((X(s,:)-center).*std_vec)/(norm(X(s,:)-center)*norm(std_vec)));
inner_points = demarcation(X(inner_nodes,:),center,(X(s,:)+X(t,:))/2,std_vec,1,demarcation_ratio,angle_bound);
outer_points = demarcation(X(outer_nodes,:),center,(X(s,:)+X(t,:))/2,std_vec,-1,demarcation_ratio,angle_bound);

inner_points = inner_nodes(inner_points>0);
outer_points = outer_nodes(outer_points>0);

demarcation_num = min(length(inner_points),length(outer_points));
inner_points_flag = zeros(1,length(inner_points));
outer_points_flag = zeros(1,length(outer_points));

for k = 1:demarcation_num
    min_angle = 9999;
    inner_point = 0;
    outer_point = 0;
    for i = 1:length(inner_points)
        if inner_points_flag(i) == 1
            continue;
        end
        for j = 1:length(outer_points)
            if outer_points_flag(j) == 1
                continue;
            end
            vec1 = X(inner_points(i),:) - center;
            vec2 = X(outer_points(j),:) - center;
            angle = acos(sum(vec1.*vec2)/(norm(vec1)*norm(vec2)));
            if angle < min_angle
                min_angle = angle;
                inner_point = i;
                outer_point = j;
            end
        end
    end
    inner_points_flag(inner_point) = 1;
    outer_points_flag(outer_point) = 1;
end

inner_points = inner_points.*int32(inner_points_flag);
inner_points = inner_points(inner_points~=0);
outer_points = outer_points.*int32(outer_points_flag);
outer_points = outer_points(outer_points~=0);

fprintf('Demarcation points generation done.\n');

out.km(out.km>-0.2) = -0.2;

[graph,graph_size] = build_graph(F',out.km,[],[],3);

teeth_num = demarcation_num + 1;
demarcation_lines = cell(1,demarcation_num);

for i = 1:demarcation_num
   [d,pre] = dijsktra(length(X),graph,graph_size,outer_points(i));
   demarcation_lines{i} = dijsktra([],pre,inner_points(i));
end

inner_boundary = cell(1,teeth_num);
outer_boundary = cell(1,teeth_num);
teeth_boundary = cell(1,teeth_num);

idx = 1;last = 1;
for i = 1:length(inner_nodes)
    if inner_nodes(i) == inner_points(idx)
        inner_boundary{idx} = inner_nodes(last:i);
        idx = idx + 1;
        last = i+1;
        if idx > demarcation_num
            break;
        end
    end
end
inner_boundary{demarcation_num+1} = inner_nodes(last:end);

idx = 1;last = 1;
for i = 1:length(outer_nodes)
    if outer_nodes(i) == outer_points(idx)
        outer_boundary{idx} = outer_nodes(last:i);
        idx = idx + 1;
        last = i + 1;
        if idx > demarcation_num
            break;
        end
    end
end
outer_boundary{demarcation_num+1} = outer_nodes(last:end);

teeth_boundary{1} = [outer_boundary{1}(end:-1:1),demarcation_lines{1},inner_boundary{1}];
for i = 2:demarcation_num
    teeth_boundary{i} = [demarcation_lines{i-1},outer_boundary{i},demarcation_lines{i}(end:-1:1),inner_boundary{i}(end:-1:1)];
end
teeth_boundary{demarcation_num+1} = [demarcation_lines{demarcation_num},inner_boundary{demarcation_num+1}(end:-1:1),outer_boundary{demarcation_num+1}];

%teeth_boundary_bak = teeth_boundary;

teeth_centers = zeros(teeth_num,3);
for i = 1:teeth_num
    teeth_centers(i,:) = sum(X(teeth_boundary{i},:))/length(teeth_boundary{i});
end

teeth_boundary_X = cell(1,teeth_num);
for i = 1:teeth_num
    teeth_boundary_X{i} = X(teeth_boundary{i},:)';
end

inner_boundary_X = cell(1,teeth_num);
for i = 1:teeth_num
    inner_boundary_X{i} = X(inner_boundary{i},:)';
end

outer_boundary_X = cell(1,teeth_num);
for i = 1:teeth_num
    outer_boundary_X{i} = X(outer_boundary{i},:)';
end

fprintf('Full contour generation done.\n');

graph = int32(graph);

teeth_centers = zeros(1,teeth_num);
teeth_centers_X = zeros(teeth_num,3);
min_dists = zeros(1,teeth_num);
for i = 1:teeth_num
    min_dists(i) = 9999;
    teeth_centers_X(i,:) = sum(X(teeth_boundary{i},:)) / length(teeth_boundary{i});
end
for i = 1:length(X)
    for j = 1:teeth_num
        temp = norm(X(i,:)-teeth_centers_X(j,:));
        if temp < min_dists(j)
            min_dists(j) = temp;
            teeth_centers(j) = i;
        end
    end
end

for i = 1:teeth_num
    while 1
        [temp_graph,temp_graph_size] = build_graph(F',out.km,double(teeth_boundary{i}),[length(teeth_boundary{i})],4);
        [d,pre] = dijsktra(length(X),temp_graph,temp_graph_size,teeth_centers(i));
        if i~=teeth_num
            that = teeth_centers(i+1);
        else
            that = teeth_centers(i-1);
        end
        path = dijsktra([],pre,that);
        if d(that) < 9999
            bias = floor(length(path)/5);
            teeth_boundary{i} = [teeth_boundary{i} path(bias:length(path)-bias)];
        else
            break;
        end
    end
end

permutation = zeros(1,teeth_num);
cnt = 1;loop = 1;
while cnt <= teeth_num
    idx = loop;
    while idx <= teeth_num
        permutation(idx) = cnt;
        cnt = cnt + 1;
        idx = idx + 3;
    end
    loop = loop + 1;
end

color = int32(zeros(1,length(X)));
for i = 1:teeth_num
    temp = dye(graph,graph_size,length(X),int32(teeth_boundary{i}),teeth_centers(i),permutation(i));
    tinted = length(temp(temp~=0));
    if tinted < length(X)/10 && tinted > length(X)/100 
        color = color + int32(color==0).*temp;
    end
    % fprintf('%d\n',tinted);
end

boundary_flag = zeros(1,length(X));
for i = 1:teeth_num
   for j = 1:length(teeth_boundary{i})
       boundary_flag(teeth_boundary{i}(j)) = 1;
   end
end

frag_color = zeros(1,length(F));

for i = 1:length(F)
    x = F(i,1);
    y = F(i,2);
    z = F(i,3);
    for j = 1:3
        if boundary_flag(F(i,j)) ~= 1
            frag_color(i) = color(F(i,j));
            break;
        end
    end
    if sum(boundary_flag(F(i,:))) == 3 && color(F(i,1)) == color(F(i,2)) && color(F(i,2)) == color(F(i,3))
        frag_color(i) = color(F(i,1));
    end
end

fprintf('Teeth annotation done.\n');

end





function points=demarcation(nodes,center1,center2,std_vec,sign,pace_ratio,angle_bound)
    len = length(nodes);
    pace = floor(len/pace_ratio);%45
    dists = zeros(1,len);
    points = zeros(1,len);
    center = center1;
    center1 = center + (center - center2);
    for i = 1:len
       vec = nodes(i,:)-center;
       angle = acos(sum(vec.*std_vec)/(norm(vec)*norm(std_vec)));
       ratio = angle / angle_bound;
       temp_center = (1-ratio)*center1 + ratio*center2;
       dists(i) = norm(nodes(i,:)-temp_center);
    end
    
    %center = (center1 + center2)/2;
    for i = pace+1:len-pace-1
        flag = 1;
        for j = i-pace:i+pace
           if (dists(i)-dists(j)) * sign < 0
                flag = -1;break;
           end
        end
        %val1 = acos(sum((nodes(i,:)-center).*(s-center))/(norm(nodes(i,:)-center)*norm(s-center)));
        %val2 = acos(sum((nodes(i,:)-center).*(t-center))/(norm(nodes(i,:)-center)*norm(t-center)));
        %angle_bias_bound = 0.2;
        if flag > 0 %&& val1 > angle_bias_bound && val2 > angle_bias_bound
            points(i) = 1;
        end
    end
    
    flag = 1;
    while flag == 1
        flag = 0;
        dems = (1:length(points)).*points;
        dems = dems(dems~=0);
        angle_diff = zeros(1,length(dems));
        for i = 2:length(dems)-1
            left = nodes(dems(i-1),:);
            right = nodes(dems(i+1),:);
            mid = nodes(dems(i),:);
            val1 = acos(sum((mid-center).*(left-center))/(norm(mid-center)*norm(left-center)));
            val2 = acos(sum((mid-center).*(right-center))/(norm(mid-center)*norm(right-center)));
            angle_diff(i) = val1 + val2;
        end
        for i = 2:length(dems)-1
            if (angle_diff(i-1)+angle_diff(i+1)-2*angle_diff(i))/angle_diff(i) > 0.5
                points(dems(i)) = 0;
                flag = 1;
                %fprintf('drop:%d\n',dems(i));
            end
        end
    end
    
    dems = (1:length(points)).*points;
    dems = dems(dems~=0);
    angle_diff = zeros(1,length(dems));
    for i = 2:length(dems)
        idx = floor((dems(i)+dems(i-1))/2);
        vec = nodes(idx,:)-center;
        angle = acos(sum(vec.*std_vec)/(norm(vec)*norm(std_vec)));
        ratio = angle / angle_bound;
        temp_center = (1-ratio)*center1 + ratio*center2;
        vec1 = nodes(dems(i-1),:) - temp_center;
        vec2 = nodes(dems(i),:) - temp_center;
        angle_diff(i) = acos(sum((vec1).*(vec2))/(norm(vec1)*norm(vec2)));
    end
    angle_average = sum(angle_diff) / (length(dems)-1);
    for i = 2:length(dems)
        if angle_diff(i) > 1.75*angle_average
            idx = floor((dems(i-1)+dems(i))/2);
            points(idx) = 1;
            %fprintf('add:%d\n',idx);
        end
    end
end

function cusps = detect_cusp(F,X,anchor_points)
    DEBUG = 1;

    if DEBUG == 1
        filepath = 'lower_right.obj';
        anchor_points = [30.7 -37.64 -1.21;-33.35 -36.14 -2.33;0.1195 7.232 -3.0];%lower_right
        [F,X,status] = readSmf(filepath);
    end

    alpha = 0.6;
    curv_upper_bound = 0;
    curv_lower_bound = -2;

    plane = [anchor_points(:,1:2) [1 1 1]'] \ anchor_points(:,3);
    plane_dist = X(:,1)*plane(1)+X(:,2)*plane(2)+plane(3)-X(:,3);

    out = tricurv_v01(F,X);
    out.km(out.km>curv_upper_bound) = curv_upper_bound;
    out.km(out.km<curv_lower_bound) = curv_lower_bound;

    H = (1-alpha)*(-out.km) + alpha*(-plane_dist);
    H_threshold = 0.625*max(H);
    H(isnan(H)) = min(H);

    [G,G_size] = build_graph((int32(F))',1+zeros(length(X),1),[],[],3);
    G = int32(G');

    is_cusp = zeros(length(X),1);
    for i = 1:length(is_cusp)
        if is_cusp(i) == -1 || H(i) < H_threshold
            is_cusp(i) = -1;
            continue;
        end
        % fprintf('%d\n',i);
        q = zeros(1000,2);
        visited = containers.Map('KeyType','int32','ValueType','int32');
        head = 1;
        tail = 1;
        q(head,:) = [i 1];
        visited(i) = 1;
        head = head + 1;
        is_valid = 0;
        while tail < head
            node = q(tail,:);
            if node(2) >= 5
                tail = tail + 1;
                is_valid = 1;
                continue;
            end
            for k = 1:G_size(node(1))
                child = G(node(1),k*2-1);
                if visited.isKey(child) == 1
                    continue;
                end
                visited(child) = 1;
                q(head,:) = [child node(2)+1];
                head = head + 1;
            end
            if H(node(1)) > H(i) || plane_dist(node(1)) < plane_dist(i)
                is_cusp(i) = -1;
                break;
            end
            if node(1) ~= i
                is_cusp(node(1)) = -1;
            end
            tail = tail + 1;
        end
        if is_valid == 0
            is_cusp(i) = -1;
        end
    end

    cusps = (1:length(is_cusp)) .* (is_cusp>-1)';
    cusps = int32(cusps(cusps~=0));

    if DEBUG == 1
        figure; 
        trisurf(F,X(:,1),X(:,2),X(:,3),H),
        colorbar vert, shading flat
        caxis([(min(H)+4*max(H))/5 max(H)])
        axis equal;
        hold on;
        scatter3(X(cusps,1),X(cusps,2),X(cusps,3),'r');
        hold off;
    end
end
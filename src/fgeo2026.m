function cmplx = fgeo2026(cmplx)
% Geometry and Forman volumes of the complex K.
cmplx(1).bfvol = ones(cmplx(1).num(1).val,1);
cmplx(1).dfvol = ones(cmplx(1).num(1).val,1);

cmplx(2).bc = zeros(cmplx(2).num(2).val,3);
cmplx(2).dir = zeros(cmplx(2).num(2).val,3);
cmplx(2).bfvol = zeros(cmplx(2).num(2).val,1);
cmplx(2).dfvol = zeros(cmplx(2).num(2).val,1);
for i=1:cmplx(2).num(2).val
    cmplx(2).bc(i,:) = mean(cmplx(1).bc(cmplx(2).bndop(1).indx(i,1:2),:));

    cmplx(2).dir(i,:) = ...
        cmplx(2).bndop(1).sgn(i,1).*...
        cmplx(1).bc(cmplx(2).bndop(1).indx(i,1),:) + ...
        cmplx(2).bndop(1).sgn(i,2).*...
        cmplx(1).bc(cmplx(2).bndop(1).indx(i,2),:);

    cmplx(2).bfvol(i) = norm(cmplx(2).dir(i,:));
    cmplx(2).dfvol(i) = norm(cmplx(2).dir(i,:));

    cmplx(2).dir(i,:) = cmplx(2).dir(i,:) ./ cmplx(2).bfvol(i);
end

cmplx(3).bc = zeros(cmplx(3).num(3).val,3);
cmplx(3).dir = zeros(cmplx(3).num(3).val,3);
cmplx(3).bfvol = zeros(cmplx(3).num(3).val,1);
cmplx(3).dfvol = zeros(2*cmplx(3).num(3).val,1);
for i=1:cmplx(3).num(3).val
    cmplx(3).bc(i,:) = mean(cmplx(1).bc(cmplx(3).bndop(1).indx(i,1:4),:));

    cnt = 0;
    for j=1:cmplx(3).num(2).val(i)
        k = cmplx(3).bndop(2).indx(i,j);

        h = norm(cross(cmplx(3).bc(i,:) - ...
            cmplx(2).bc(k,:),cmplx(2).dir(k,:)));
        cmplx(3).bfvol(i) = cmplx(3).bfvol(i) + 0.5*cmplx(2).bfvol(k)*h;

        if j==1 || j==4
            cmplx(3).dfvol(2*(i-1)+1) = ...
                cmplx(3).dfvol(2*(i-1)+1) + 0.5*cmplx(2).dfvol(k)*h;
        else
            cmplx(3).dfvol(2*(i-1)+2) = ...
                cmplx(3).dfvol(2*(i-1)+2) + 0.5*cmplx(2).dfvol(k)*h;
        end

        v1 = cmplx(3).bndop(2).sgn(i,j)*cmplx(2).dir(k,:);
        vv1 = cmplx(3).bndop(2).sgn(i,j)*cmplx(2).bndop(1).indx(k,:).*cmplx(2).bndop(1).sgn(k,:);
        vv1e = vv1(vv1>0);
        for j2=1:cmplx(3).num(2).val(i)
            k2 = cmplx(3).bndop(2).indx(i,j2);
            if k==k2
                continue
            end
            v2 = cmplx(3).bndop(2).sgn(i,j2)*cmplx(2).dir(k2,:);
            vv2 = cmplx(3).bndop(2).sgn(i,j2)*cmplx(2).bndop(1).indx(k2,:).*cmplx(2).bndop(1).sgn(k2,:);
            if any(vv2==-vv1e)
                cnt = cnt + 1;
                tmp = cross(v1,v2)./norm(cross(v1,v2));
                cmplx(3).dir(i,:) = cmplx(3).dir(i,:) + tmp;
                break
            elseif any(vv2==vv1e)
                break
            end
        end

    end
    cmplx(3).dir(i,:) = cmplx(3).dir(i,:) / cnt;
    cmplx(3).dir(i,:) = cmplx(3).dir(i,:)./norm(cmplx(3).dir(i,:));

end

cmplx(4).bfvol = zeros(cmplx(4).num(4).val,1);
for i=1:cmplx(4).num(4).val
    cmplx(4).bc(i,:) = mean(cmplx(1).bc(cmplx(4).bndop(1).indx(i,:),:));
    for j=1:cmplx(4).num(3).val(i)
        k = cmplx(4).bndop(3).indx(i,j);
        for l=1:cmplx(3).num(2).val(k)
            m=cmplx(3).bndop(2).indx(k,l);
            edgs=[cmplx(3).bc(k,:);...
                cmplx(1).bc(cmplx(2).bndop(1).indx(m,1:2),:)] - ...
                cmplx(4).bc(i,:);
            cmplx(4).bfvol(i) = cmplx(4).bfvol(i) + ...
                cmplx(4).bndop(3).sgn(i,j)*cmplx(3).bndop(2).sgn(k,l)/6*...
                (dot(edgs(1,:),cmplx(2).bndop(1).sgn(m,1)*...
                cross(edgs(2,:),edgs(3,:))));
        end
    end
end
cmplx(4).bfvol=abs(cmplx(4).bfvol);

cmplx(1).fvol  = zeros(cmplx(1).num(1).val,1);
for i=1:cmplx(1).num(1).val
    for j=1:cmplx(1).num(4).val(i)
        k = cmplx(1).bndop(4).indx(i,j);

        cmplx(1).fvol(i)  = cmplx(1).fvol(i)  + cmplx(4).bfvol(k);
    end
    cmplx(1).fvol(i)  = cmplx(1).fvol(i)  / (2^3*cmplx(1).bfvol(i));
end

cmplx(2).fvol  = zeros(cmplx(2).num(2).val,1);
for i=1:cmplx(2).num(2).val
    cnt = 0;
    found = [];
    e_vrts = nonzeros(cmplx(2).bndop(1).indx(i,:));
    for j=nonzeros(cmplx(2).bndop(4).indx(i,:))'
        for k = nonzeros(cmplx(4).bndop(3).indx(j,:))'
            f_vrts = nonzeros(cmplx(3).bndop(1).indx(k,:));
            if sum(e_vrts==f_vrts','all')==1
                if ~any(found==k)
                    cmplx(2).fvol(i) = cmplx(2).fvol(i) + cmplx(3).bfvol(k);
                    found = [found,k];
                    cnt = cnt+1;
                end
            end
        end
    end
    cmplx(2).fvol(i)  = cmplx(2).fvol(i)  / (2^3*cmplx(2).bfvol(i));
end

cmplx(3).fvol  = zeros(cmplx(3).num(3).val,1);
for i=1:cmplx(3).num(3).val
    found = [];
    f_vrts = nonzeros(cmplx(3).bndop(1).indx(i,:));
    for j=nonzeros(cmplx(3).bndop(4).indx(i,:))'
        for k = nonzeros(cmplx(4).bndop(2).indx(j,:))'
            e_vrts = nonzeros(cmplx(2).bndop(1).indx(k,:));
            if sum(e_vrts==f_vrts','all')==1
                if ~any(found==k)
                    cmplx(3).fvol(i) = cmplx(3).fvol(i) + cmplx(2).bfvol(k);
                    found = [found,k];
                    cnt = cnt+1;
                end
            end
        end
    end
    cmplx(3).fvol(i)  = cmplx(3).fvol(i)  / (2^3*cmplx(3).bfvol(i));
end

cmplx(4).fvol  = zeros(cmplx(4).num(4).val,1);
for i=1:cmplx(4).num(4).val
    for j=1:cmplx(4).num(1).val(i)
        k = cmplx(4).bndop(1).indx(i,j);

        cmplx(4).fvol(i) = cmplx(4).fvol(i) + cmplx(1).bfvol(k);
    end
    cmplx(4).fvol(i)  = cmplx(4).fvol(i)  / (2^3*cmplx(4).bfvol(i));
end

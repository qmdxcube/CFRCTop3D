function plotOptimizedDesign(hx, hy, hz, xPhys, yPhys, pxPhys, pyPhys, pzPhys, U)

[nelz, nely, nelx] = size(xPhys);
nNode = (nelx+1)*(nely+1)*(nelz+1);

[Z, Y, X] = ndgrid(0:nelz, 0:nely, 0:nelx);
coords = [X(:)*hx, Y(:)*hy, Z(:)*hz];

nodenrs = int32(reshape(1:nNode,1+nelz,1+nely,1+nelx));
edofVec = reshape(nodenrs(1:end-1, 1:end-1, 1:end-1), [], 1);
nids = edofVec+int32([0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1,1+[0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1]]);
facePattern = [1 4 8 5; 4 8 7 3; 3 7 6 2; 2 6 5 1; 1 4 3 2; 8 5 6 7];

xthreshold = 0.5; ythreshold = 0.01;
validElems = (xPhys >= xthreshold) & (yPhys > ythreshold);
if ~any(validElems(:)), warning('No elements satisfy the threshold'); return; end

validNids = nids(validElems(:), :);
elementFace = reshape(validNids(:, facePattern')', 4, [])';

for k=1:size(U,2)
    % U_mag
    Ux = U(1:3:end,k); Uy = U(2:3:end,k); Uz = U(3:3:end,k);
    U_show = sqrt(Ux.^2 + Uy.^2 + Uz.^2);

    % ---Displacement Mode---
    figure('Color', 'w', 'Name', 'Displacement Field');
    patch('Faces', elementFace, 'Vertices', coords, ...
        'FaceVertexCData', U_show, 'FaceColor', 'interp', ...
        'EdgeAlpha', 0.1, 'EdgeColor', 'k');
    setup_3d_style('Displacement Magnitude', [min(U_show), max(U_show)], jet);
end

% --- Fiber Design & yPhys ---
figure('Color', 'w', 'Name', 'Fiber Design');
set(gcf, 'Renderer', 'opengl'); 
hold on;

patch('Faces', elementFace, 'Vertices', coords, ...
    'FaceColor', [0.85, 0.85, 0.9], 'FaceAlpha', 0.2, ...
    'EdgeAlpha', 0.05, 'EdgeColor', 'k');

[CZ, CY, CX] = ndgrid((0.5:nelz-0.5)*hz, (0.5:nely-0.5)*hy, (0.5:nelx-0.5)*hx);
scale = hx * 0.8;
px = pxPhys(validElems) * scale;
py = pyPhys(validElems) * scale;
pz = pzPhys(validElems) * scale;

q = quiver3(CX(validElems)-px/2, CY(validElems)-py/2, CZ(validElems)-pz/2, px, py, pz, ...
    'AutoScale', 'off', 'LineWidth', 1.2, 'ShowArrowHead', 'off');

cmap_fiber = jet(256);
y_min = min(yPhys(:)); y_max = max(yPhys(:));
if y_max == y_min, y_max = y_min + 1e-6; end
norm_y = (yPhys(validElems) - y_min) / (y_max - y_min);
ind = floor(norm_y * 255) + 1;
fiber_colors = uint8([repelem(cmap_fiber(ind,:), 2, 1), ones(numel(ind)*2, 1)].' * 255);
set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', fiber_colors);
setup_3d_style('Fiber Volume Fraction', [y_min, y_max], jet);
hold off;
end
function MESH = generateSquareMesh(xL,xR,yL,yR,IMAX,JMAX,periodic)

global R
%GENERATESQUAREMESH Summary of this function goes here
%   Detailed explanation goes here
MESH = struct();
%MESH.nNode = (IMAX+1)*(JMAX+1); % number of nodes
%MESH.nElem = 2*IMAX*JMAX;

%MESH.x = zeros(MESH.nNode,2);

dx = (xR-xL)/IMAX;
dy = (yR-yL)/JMAX;

count = 0;
% internal nodes have random displacements
for j = 1:JMAX+1
    %if (mod(j,2) == 0)
    %    count = count + 1;
    %    MESH.x(count,1) = xL;
    %    MESH.x(count,2) = yL+ (j-1)*dy;
    %    for i = 1:IMAX
    %        count = count + 1;
    %        MESH.x(count,1) = xL+ 0.5*dx + (i-1)*dx;
    %        MESH.x(count,2) = yL+ (j-1)*dy;
    %    end
    %    count = count + 1;
    %    MESH.x(count,1) = xR;
    %    MESH.x(count,2) = yL+ (j-1)*dy;
    %else
        for i = 1:IMAX + 1
            count = count + 1;
            MESH.x(count,1) = xL + (i-1)*dx;
            MESH.x(count,2) = yL+ (j-1)*dy;
        end
   % end



end

% rotate mesh 
%center = 0.5*[xL+xR; yL + yR];
%xNew = zeros(count,2); 
%for j = 1:count 
%    aux = center + R*( MESH.x(j,:)' - center);
%    xNew(j,:) = aux';
%end
%MESH.x = xNew;

%count = count +1;
%x1 = xR - dx;
%x2 = xR;
%x3 = xR;
%y1 = yR;
%y2 = yR;
%y3 = yR-dy;
%MESH.x(count,1) = 1/2*(x1 + x2);
%MESH.x(count,2) = 1/2*(y2 + y3);
%count = count + 1;
%MESH.x(count,1) = 1/2*(x1 + x2) -dx;
%MESH.x(count,2) = 1/2*(y2 + y3);
MESH.nNode = count;

MESH.tri = delaunay(MESH.x(:,1),MESH.x(:,2));
MESH.nElem = size(MESH.tri,1);

MESH.verini2verper = 1:MESH.nNode;
% we simply change the connectivity matrix tri and fuse some nodes of the
% grid
%
% generate the list of nodes that need to be removed
MESH.nRemoveX = 0;
MESH.nRemoveY = 0;
MESH.nFuseX   = 0;
MESH.nFuseY   = 0;
if(periodic(1)==1)
    % if we are periodic in x, remove all nodes on the right boundary and
    % fuse them into the nodes on the left boundary
    for iNode=1:MESH.nNode
        if( abs(MESH.x(iNode,1)-xR)<1e-13 )
            MESH.nRemoveX = MESH.nRemoveX + 1;
            MESH.RemoveX(MESH.nRemoveX) = iNode;
        end
        if( abs(MESH.x(iNode,1)-xL)<1e-13 )
            MESH.nFuseX = MESH.nFuseX + 1;
            MESH.FuseX(MESH.nFuseX) = iNode;
        end
    end
    [ys,id]=sort(MESH.x(MESH.FuseX,2));
    MESH.FuseX=MESH.FuseX(id);
    [ys,id]=sort(MESH.x(MESH.RemoveX,2));
    MESH.RemoveX=MESH.RemoveX(id);
    if(MESH.nFuseX ~= MESH.nRemoveX)
        disp(sprintf('number of x nodes incompatible %d (L) and %d (R)', MESH.nFuseX, MESH.nRemoveX));
        return
    end
end
if(periodic(2)==1)
    % if we are periodic in y, remove all nodes on the top boundary and
    % fuse them into the nodes on the bottom boundary
    for iNode=1:MESH.nNode
        if( abs(MESH.x(iNode,2)-yR)<1e-13 )
            MESH.nRemoveY = MESH.nRemoveY + 1;
            MESH.RemoveY(MESH.nRemoveY) = iNode;
        end
        if( abs(MESH.x(iNode,2)-yL)<1e-13 )
            MESH.nFuseY = MESH.nFuseY + 1;
            MESH.FuseY(MESH.nFuseY) = iNode;
        end
    end
    [ys,id]=sort(MESH.x(MESH.FuseY,1));
    MESH.FuseY=MESH.FuseY(id);
    [ys,id]=sort(MESH.x(MESH.RemoveY,1));
    MESH.RemoveY=MESH.RemoveY(id);
    if(MESH.nFuseY ~= MESH.nRemoveY)
        disp(sprintf('number of y nodes incompatible %d (L) and %d (R)', MESH.nFuseY, MESH.nRemoveY));
        return
    end
end
%
nodemap = 1:MESH.nNode;
tri = MESH.tri;
if(periodic(1)==1)
    % if we have periodic boundary conditions in x direction
    for j=1:MESH.nRemoveX
        iRemoveNode = MESH.RemoveX(j);     % global number of the node that must be removed
        iFuseNode   = MESH.FuseX(j);       % global number of the node to which it is fused
        for iElem=1:MESH.nElem
            for kLocNode=1:3 % run over local nodes
                iNode = tri(iElem,kLocNode);
                if(iNode==iRemoveNode)
                    tri(iElem,kLocNode) = iFuseNode;
                    nodemap(iRemoveNode) = iFuseNode;
                end
            end
        end
    end
end


if(periodic(2)==1)
    % if we have periodic boundary conditions in x direction
    for i=1:MESH.nRemoveY
        iRemoveNode = nodemap(MESH.RemoveY(i));      % global number of the node that must be removed
        iFuseNode   = nodemap(MESH.FuseY(i));        % global number of the node to which it is fused
        for iElem=1:MESH.nElem
            for kLocNode=1:3
                iNode = tri(iElem,kLocNode);
                if( i<MESH.nRemoveY && MESH.tri(iElem,kLocNode)==MESH.RemoveY(MESH.nRemoveY))
                    continue
                end
                if(iNode==iRemoveNode)
                    tri(iElem,kLocNode) = iFuseNode;
                    nodemap(MESH.RemoveY(i)) = iFuseNode;
                end
            end
        end
    end
end

MESH.verini2verper = nodemap;


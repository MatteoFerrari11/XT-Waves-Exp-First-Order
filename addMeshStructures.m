function MESH = addMeshStructures(MESH)

%global R

R = [0 1; -1 0];

% compute volume of elements
MESH.vol = zeros(MESH.nElem,1);
MESH.xb = zeros(MESH.nElem,2);
MESH.T = zeros(MESH.nElem,2,2);
MESH.coord = zeros(MESH.nElem,2,3);
%
for iElem = 1:MESH.nElem
    coord = MESH.x(MESH.tri(iElem,:),:)';
    MESH.coord(iElem,:,:) = coord;
    T = [coord(:,2) - coord(:,1), coord(:,3) - coord(:,1)];
    MESH.Tmat(iElem,:,:) = T;
    MESH.vol(iElem) = 0.5 * det(T);
    MESH.xb(iElem,:) = 1./3 * sum(coord,2);
end

% start with voronoi dual
MESH.nNode2Element = zeros(MESH.nNode,1);
MESH.node2Element = zeros(MESH.nNode,10);
MESH.node2Element2locNode = zeros(MESH.nNode,10);

for iElem = 1:MESH.nElem
    for iLocNode = 1:3
        iNode = MESH.verini2verper(MESH.tri(iElem,iLocNode));
        MESH.nNode2Element(iNode) = MESH.nNode2Element(iNode) + 1;
        MESH.node2Element(iNode,MESH.nNode2Element(iNode)) = iElem;
        MESH.node2Element2locNode(iNode,MESH.nNode2Element(iNode)) = iLocNode;
    end
end

maxnface = 3*MESH.nElem;
MESH.FaceDef = [2 3;
    3 1;
    1 2];
MESH.Neighbor = zeros(MESH.nElem,3);
MESH.NeighborLocFace = zeros(MESH.nElem,3);
MESH.nFace = 0;
%
MESH.triFace = zeros(MESH.nElem,3);
MESH.face2Element = zeros(maxnface,2);
MESH.face2ElementLocFace = zeros(maxnface,2);
MESH.faces = zeros(maxnface,2);
MESH.element2faceLocNodes = zeros(MESH.nElem, 3, 2);
MESH.element2faceSign = zeros(MESH.nElem,3);
% look for faces and neighbors
for iElem = 1:MESH.nElem
    for iLocFace = 1:3
        iNode1 = MESH.verini2verper(MESH.tri(iElem,MESH.FaceDef(iLocFace,1)));
        iNode2 = MESH.verini2verper(MESH.tri(iElem,MESH.FaceDef(iLocFace,2)));
        %
        if (MESH.triFace(iElem,iLocFace) == 0) % if the face has not been added yet
            MESH.nFace = MESH.nFace + 1;
            MESH.faces(MESH.nFace,:) = [iNode1,iNode2];
            MESH.triFace(iElem,iLocFace) = MESH.nFace;
            MESH.face2Element(MESH.nFace,1) = iElem;
            MESH.face2ElementLocFace(MESH.nFace,1) = iLocFace;

            MESH.element2faceLocNodes(iElem,iLocFace,1) = MESH.FaceDef(iLocFace,1);
            MESH.element2faceLocNodes(iElem,iLocFace,2) = MESH.FaceDef(iLocFace,2);
            MESH.element2faceSign(iElem,iLocFace) = 1;
            %Mi
            for jDual = 1:MESH.nNode2Element(iNode1)
                jElem = MESH.node2Element(iNode1,jDual);
                for jLocFace = 1:3
                    jNode1 = MESH.verini2verper(MESH.tri(jElem,MESH.FaceDef(jLocFace,1)));
                    jNode2 = MESH.verini2verper(MESH.tri(jElem,MESH.FaceDef(jLocFace,2)));
                    if (iNode1 == jNode2 && iNode2 == jNode1)
                        MESH.Neighbor(iElem,iLocFace) = jElem;
                        MESH.Neighbor(jElem,jLocFace) = iElem;
                        MESH.NeighborLocFace(iElem,iLocFace) = jLocFace;
                        MESH.NeighborLocFace(jElem,jLocFace) = iLocFace;
                        %
                        MESH.triFace(jElem,jLocFace) = MESH.nFace;
                        MESH.face2Element(MESH.nFace,2) = jElem;
                        MESH.face2ElementLocFace(MESH.nFace,2) = jLocFace;
                        %
                        MESH.element2faceLocNodes(jElem,jLocFace,1) = MESH.FaceDef(jLocFace,2);
                        MESH.element2faceLocNodes(jElem,jLocFace,2) = MESH.FaceDef(jLocFace,1);
                        MESH.element2faceSign(jElem,jLocFace) = -1;
                    end
                end
            end
        end
    end
end
% set the right dimension
MESH.faces = MESH.faces(1:MESH.nFace,:);
MESH.face2Element = MESH.face2Element(1:MESH.nFace,:);
MESH.face2ElementLocFace = MESH.face2ElementLocFace(1:MESH.nFace,:);


% node2face
MESH.nNode2face = zeros(MESH.nNode);
MESH.node2face = zeros(MESH.nNode, 10);
for iFace = 1:MESH.nFace
    for iLocNode = 1:2 
        iNode = MESH.faces(iFace,iLocNode);
        MESH.nNode2face(iNode) = MESH.nNode2face(iNode) + 1;
        MESH.node2face(iNode, MESH.nNode2face(iNode)) = iFace;
    end
end

% compute normal vectors and barycenter of faces

MESH.element2FaceXb = zeros(MESH.nElem,3,2);
MESH.element2FaceNormal = zeros(MESH.nElem,3,2); % outer normal
MESH.faceNormal = zeros(MESH.nFace,2); % face oriented normal
MESH.faceXb = zeros(MESH.nFace,2);
MESH.faceTangent = zeros(MESH.nFace,2);

for iElem = 1:MESH.nElem 
    for iLocFace = 1:3 
        iFace = MESH.triFace(iElem,iLocFace);

        locNode1 = MESH.FaceDef(iLocFace,1);
        locNode2 = MESH.FaceDef(iLocFace,2);

        x1 = MESH.x(MESH.tri(iElem,locNode1),:);
        x2 = MESH.x(MESH.tri(iElem,locNode2),:);
        %
        MESH.element2FaceXb(iElem,iLocFace,:) = 0.5*(x1(:) + x2(:));
        MESH.faceXb(iFace,:) = 0.5*(x1(:) + x2(:));
        % 
        t = x2(:) - x1(:);
        eta = R*t;
        MESH.element2FaceNormal(iElem,iLocFace,:) = eta(:);
        if (MESH.element2faceSign(iElem,iLocFace) == 1)
            MESH.faceTangent(iFace,:) = t(:);
            MESH.faceNormal(iFace,:) = eta(:);
        else
            MESH.faceNormal(iFace,:) = -eta(:);
            MESH.faceTangent(iFace,:) = -t(:);
        end
    end
end

% compute boundary conditions for the LDC 
MESH.bndFlagFace = zeros(MESH.nFace,1);
MESH.bndFlagNode = zeros(MESH.nNode,1);
MESH.global2InternalFace = zeros(MESH.nFace,1);
MESH.internal2GlobalFace = [];
count = 0;
for iFace = 1:MESH.nFace
    %xb(:) = MESH.faceXb(iFace,:);
    %if (abs(xb(1) - xL)<1e-13 || abs(xb(1) - xR)<1e-13  || abs(xb(2) - yL)<1e-13 || abs(xb(2) - yR)<1e-13 )
    %    MESH.bndFlagFace(iFace) = 1;
    %end
    rElem = MESH.face2Element(iFace,2); 
    if (rElem == 0) 
        MESH.bndFlagFace(iFace) = 1;
        for ii = 1:2
            iNode = MESH.faces(iFace,ii);
            MESH.bndFlagNode(iNode) = 1;
        end
    else
        count = count + 1;
        MESH.internal2GlobalFace(count) = iFace; 
        MESH.global2InternalFace(iFace) = count;
    end
end
MESH.nInternalFace = count;



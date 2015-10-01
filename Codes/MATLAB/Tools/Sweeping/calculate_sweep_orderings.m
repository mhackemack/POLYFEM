%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate Sweep Chunks
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = calculate_sweep_orderings(data, mesh)
global glob
t = tic;
if glob.print_info
    disp('-> Begin Sweep Ordering Computations.')
end
% Process some input information
% ------------------------------
data = calculate_ND_sweep_orderings(data, mesh);
if glob.print_info
    disp(['-> Total Sweep Ordering Time:  ',num2str(toc(t))])
    disp(' ')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = calculate_1D_sweep_orderings(data, mesh)
sweep.NumberSweepChunks = data.Neutronics.Transport.AngleSets;
% Loop through angle sets
for m=1:length(data.Neutronics.Transport.AngleSets)
    ave_ang = data.Neutronics.Transport.AverageAngles(m,:);
    if ave_ang > 0
        sweep.CellSweepOrder{m} = 1:mesh.TotalCells;
        sweep.DownstreamCells{m} = cell(mesh.TotalCells, 1);
        sweep.DownstreamFaces{m} = cell(mesh.TotalCells, 1);
        sweep.UpstreamCells{m} = cell(mesh.TotalCells, 1);
        sweep.UpstreamFaces{m} = cell(mesh.TotalCells, 1);
        sweep.CycleFree(m) = true;
        for c=1:mesh.TotalCells
            if c == 1
                sweep.DownstreamCells{m}{c} = c+1;
                sweep.UpstreamCells{m}{c} = 1;
            elseif c == mesh.TotalCells
                sweep.DownstreamCells{m}{c} = mesh.TotalCells;
                sweep.UpstreamCells{m}{c} = c-1;
            else
                sweep.DownstreamCells{m}{c} = c+1;
                sweep.UpstreamCells{m}{c} = c-1;
            end
            sweep.DownstreamFaces{m}{c} = c+1;
            sweep.UpstreamFaces{m}{c} = c;
        end
    else
        sweep.CellSweepOrder{m} = mesh.TotalCells:-1:1;
        sweep.DownstreamCells{m} = cell(mesh.TotalCells, 1);
        sweep.DownstreamFaces{m} = cell(mesh.TotalCells, 1);
        sweep.UpstreamCells{m} = cell(mesh.TotalCells, 1);
        sweep.UpstreamFaces{m} = cell(mesh.TotalCells, 1);
        sweep.CycleFree(m) = true;
        for c=1:mesh.TotalCells
            if c == 1
                sweep.DownstreamCells{m}{c} = 1;
                sweep.UpstreamCells{m}{c} = c+1;
            elseif c==mesh.TotalCells
                sweep.DownstreamCells{m}{c} = c-1;
                sweep.UpstreamCells{m}{c} = mesh.TotalCells;
            else
                sweep.DownstreamCells{m}{c} = c-1;
                sweep.UpstreamCells{m}{c} = c+1;
            end
            sweep.DownstreamFaces{m}{c} = c;
            sweep.UpstreamFaces{m}{c} = c+1;
        end
    end
end
data.Neutronics.Transport.Sweeping = sweep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = calculate_ND_sweep_orderings(data, mesh)
sweep.NumberSweepChunks = data.Neutronics.Transport.AngleSets;
% Loop through angle sets
for m=1:length(data.Neutronics.Transport.AngleSets)
    ave_ang = data.Neutronics.Transport.AverageAngles(m,:)';
    i = []; j = []; ind = []; faces = [];
    ui = []; uj = []; uind = []; ufaces = [];
    ubf = cell(mesh.TotalCells, 1); dbf = cell(mesh.TotalCells, 1);
    for c=1:mesh.TotalCells
        cf = mesh.CellFaces{c};
        for ff=1:length(cf)
            f = cf(ff);
            fid = mesh.FaceID(f);
            fnorm = mesh.FaceNormal(f,:);
            fcells = mesh.FaceCells(f,:);
            % Handle Boundary Faces Separately
            if fid~=0
                if fnorm*ave_ang > 0
                    dbf{c} = [dbf{c},f];
                else
                    ubf{c} = [ubf{c},f];
                end
                continue;
            end
            if c ~= fcells(1)
                fcells = fcells(2:-1:1);
                fnorm = -fnorm;
            end
            if fnorm*ave_ang > 0
                i = [i;c]; j = [j;fcells(2)]; 
                ind = [ind;1]; faces = [faces;f];
            else
                ui = [ui;c]; uj = [uj;fcells(2)]; 
                uind = [uind;1]; ufaces = [ufaces;f];
            end
        end
    end
    tmc = double(mesh.TotalCells);
    i = double(i); j = double(j);
    ui = double(ui); uj = double(uj);
    mat = sparse(i, j, ind, tmc, tmc);
    umat = sparse(ui, uj, uind, tmc, tmc);
    facemat = sparse(i, j, faces, tmc, tmc);
    ufacemat = sparse(ui, uj, ufaces, tmc, tmc);
    if ~graphisdag(mat)
        bg = biograph(mat);
        if mesh.Dimension == 2
            plot_mesh(mesh,0,0);
            view(bg);
        end
%         asp = bg.allshortestpaths; asp(asp == inf) = 0;
%         sweep.ShortestPath(m) = max(max(asp));
%         sweep.CycleFree(m) = bg.isdag;
%         sweep.CellSweepOrder{m} = bg.topoorder;
        error('Cannot currently break cycles.')
    else
        sweep.CellSweepOrder{m} = topological_order(mat);
        sweep.CycleFree(m) = true;
    end
    % Calculate Downstream/UpstreamCells
    sweep.DownstreamCells{m} = cell(mesh.TotalCells, 1);
    sweep.UpstreamCells{m} = cell(mesh.TotalCells, 1);
    sweep.DownstreamFaces{m} = cell(mesh.TotalCells, 1);
    sweep.UpstreamFaces{m} = cell(mesh.TotalCells, 1);
    for c=1:mesh.TotalCells
        % downstream interior
        tcmat = mat(c,:);
        [~,cID] = find(tcmat);
        sweep.DownstreamCells{m}{c} = cID;
        sweep.DownstreamFaces{m}{c} = full(facemat(c,cID));
        % downstream boundary
        ndbf = length(dbf{c});
        sweep.DownstreamCells{m}{c} = [sweep.DownstreamCells{m}{c}, ones(1,ndbf)*c];
        sweep.DownstreamFaces{m}{c} = [sweep.DownstreamFaces{m}{c}, dbf{c}];
        % upstream interior
        tcmat = umat(c,:);
        [~,cID] = find(tcmat);
        sweep.UpstreamCells{m}{c} = cID;
        sweep.UpstreamFaces{m}{c} = full(ufacemat(c,cID));
        % upstream boundary
        nubf = length(ubf{c});
        sweep.UpstreamCells{m}{c} = [sweep.UpstreamCells{m}{c}, ones(1,nubf)*c];
        sweep.UpstreamFaces{m}{c} = [sweep.UpstreamFaces{m}{c}, ubf{c}];
    end
end
data.Neutronics.Transport.Sweeping = sweep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = findcycles(G)
numNodes = size(G,1); 
for n = 1:numNodes
   [D,P]=graphtraverse(G,n);
   for d = D
       if G(d,n)
           graphpred2path(P,d)
       end
   end
   G(n,:)=0; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
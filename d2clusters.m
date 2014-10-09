function [clusters] = d2clusters( db, k)
  %%  
  % INPUT
  % k: number of clusters
  % db.stride: size of supports
  % db.w: prob of supports
  % db.supp: supports

  % OUTPUT
  % clusters: k convergent d2 clusters
  %

  global stdoutput IDX;
  
  if nargin == 1
    k = 5;
  end

  clusters = [];
  n = length(db{1}.stride); % size of total samples
  labels = randi(k,1,n);
  

  nROUND = 10;
  centroid_init = randi(n,[k,1]);
  nphase = length(db);
  clusters = cell(k);
  for j=1:k
      for i=1:nphase          
          tmps = sum(db{i}.stride(1:centroid_init(j)-1));
          strips = tmps:tmps+db{i}.stride(centroid_init(j))-1;
          clusters{j}{i}.supp = db{i}.supp(:,strips);
          clusters{j}{i}.w = db{i}.w(strips);
      end
  end

  for i=1:nROUND+1
    fprintf(stdoutput, 'Round %d ... ', i);
    % relabel based on distance
    D=zeros(k,n,nphase);
    for p=1:nphase
        for j=1:k
            strip = 1;
            for idx = 1:n
                D(j,idx,nphase) = kantorovich(clusters{j}{p}.supp, clusters{j}{p}.w,...
                    db{p}.supp(:,strip:strip+db{p}.stride(idx)-1), ...
                    db{p}.w(strip:strip+db{p}.stride(idx)-1));
                strip = strip + db{p}.stride(idx);
            end
        end
    end
    
    coeff=ones(nphase,1);
    DC = zeros(k,n);
    for p=1:nphase
        DC = DC + coeff(p)* D(:,:,p);
    end
    labelspast = labels;
    [~, labels] = min(DC);
    fprintf(stdoutput, '%d labels change \n',sum(labelspast ~= labels));
    if (i==nROUND+1)||(sum(labelspast ~= labels) == 0) 
        % export rank to each cluster centroid
        [~, IDX] = sort(DC,2);
        break;
    end
    

    % compute the centroids
    for j=1:k
      fprintf(stdoutput, '\n\t cluster %d - ', j);
      clusters{j} = centroid(j, labels, db, clusters{j});
    end

  end
end



function [c] = centroid( lb, labels, db, c0)
% INPUT
% lb: to compute the multiphase centroid of points with label = lb

% OUTPUT
% c: multiphase centroid of selected points
  


  global stdoutput ctime bufferc;   
  nphase = length(db);
  c = cell(nphase,1);
  
  if nargin == 3
      c0 = cell(nphase,1);
  end
  
  for i=1:nphase
    fprintf(stdoutput, '\n\t\tphase %d: ', i);
    warmlabels = getwarm(labels, db{i}.stride);
    w = db{i}.w(lb == warmlabels);
    supp = db{i}.supp(:,lb==warmlabels);    
    stride = db{i}.stride(lb == labels);
    
    ctimer = tic;c{i} = centroid_sphADMM(stride, supp, w, c0{i});ctime(1)=toc(ctimer); 
    %bufferc{1} = c{i};
    %ctimer = tic;c{i} = centroid_sphLP(stride{i}, supp{i}, w{i});ctime(2)=toc(ctimer);
    %bufferc{2} = c{i};
  end
  
end


function [warmlabels] = getwarm(lbs, stride)
% INPUT
% lbs: [1,2,1,3,2]
% stride: [4,3,2,3,4]

% OUTPUT
% warmlabels: [1 1 1 1 2 2 2 1 1 3 3 3 2 2 2 2]
%

  len = sum(stride);
  warmlabels = zeros(1,len);
  pos=1;
  for j=1:length(stride)
    warmlabels(pos:pos+stride(j)-1) = lbs(j);
    pos = pos + stride(j);
  end

end




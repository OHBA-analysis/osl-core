function plot_network(self,cmat,threshold,mapping)
	% Plot connection matrices as network with edges
    if nargin < 4 || isempty(mapping)
       	mapping = @(w) (4*w).^2; % Apply a transformation to map connection weight to line width
    end
            
	if nargin < 3 || isempty(threshold) 
		threshold = 0.95;
	end
	
	if size(cmat,1) ~= size(cmat,2)
		error('Input matrix must be square');
	end

	if size(cmat,1) ~= self.n_parcels
		error(sprintf('Correlation matrix has %d ROIs, but parcellation has %d parcels',size(cmat,1),self.n_parcels));
	end
	
	roi_centers = zeros(self.n_parcels,3);
    
    [coords,weights] = self.roi_coords;
	for j = 1:self.n_parcels
		roi_centers(j,:) = sum(bsxfun(@times,coords{j},weights{j}))./sum(weights{j});
	end

	f=figure;
	mni_coords = self.template_coordinates;
	shp = alphaShape(mni_coords(:,1),mni_coords(:,2),mni_coords(:,3));
	h_brain = shp.plot;
	h_brain.FaceAlpha = 0.1;
	h_brain.FaceColor = [0.5 0.5 1];
	h_brain.EdgeColor = 'w';
	hold on

	scatter3(roi_centers(:,1),roi_centers(:,2),roi_centers(:,3),30,'bo','filled');

	cmat = abs(cmat);
	threshold = prctile(cmat(:),100*threshold);
	ind = find( triu(ones(size(cmat)),1) & cmat >= threshold);
	[from,to] = ind2sub(size(cmat),ind);

	h = [];

	for j = 1:length(from)
		c = roi_centers([from(j) to(j)],:);
		h(j) = plot3(c(:,1),c(:,2),c(:,3),'r','LineWidth',mapping(cmat(from(j),to(j))));
	end

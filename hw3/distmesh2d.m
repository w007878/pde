function [p, t] = distmesh2d(fd, fh, h0, bbox, pfix, varargin)
    %% Init
    Fnscale = 1.2; 
    deltat = .2; 
    dnptol = .001; 
    ttol = .1; 
    Geps = .00001 * h0; 
    Deps = sqrt(eps) * h0;

    %% Initial equilateral triangles
    [x,y] = meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
    x(2:2:end,:) = x(2:2:end,:) + h0/2;
    p = [x(:), y(:)];     
    p = p(feval(fd,p,varargin{:})<Geps,:); 
    r0 = 1./feval(fh,p,varargin{:}).^2; % Probability to keep point
    p = [pfix; p(rand(size(p,1),1)<r0./max(r0),:)]; % Rejection method
    N = size(p,1);
    
    %% Iteration
    pold = inf; 
    while 1
        % Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((p-pold).^2,2))/h0) > ttol 
            pold = p;
            t = delaunayn(p); % List of triangles
            pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3; % Compute centroids
            t = t(feval(fd,pmid,varargin{:})<-Geps,:); % Keep interior triangles

            bars = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])]; % Interior bars duplicated
            bars = unique(sort(bars,2),'rows'); % Bars as node pairs
        % plot
            trimesh(t, p(:,1), p(:,2), zeros(N,1))
            view(2), axis equal, axis off, drawnow
        end
        
        barvec = p(bars(:,1),:) - p(bars(:,2),:); % List of bar vectors
        L = sqrt(sum(barvec.^2,2)); % L = Bar lengths
        hbars = feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2, varargin{:});
        L0 = hbars*Fnscale*sqrt(sum(L.^2)/sum(hbars.^2)); % L0 = Desired lengths
        F = max(L0-L,0); % Bar forces (scalars)
        Fvec = F./L*[1,1].*barvec; % Bar forces (x,y components)
        Ftot = full(sparse(bars(:, [1,1,2,2]), ones(size(F))*[1,2,1,2], [Fvec,-Fvec], N, 2));
        Ftot(1:size(pfix,1), :)=0; % Force = 0 at fixed points
        p = p + deltat*Ftot; % Update node positions
        % 7. Bring outside points back to the boundary
        d = feval(fd,p,varargin{:}); 
        ix = d>0; % Find points outside (d>0)
        dgradx = (feval(fd,[p(ix,1)+Deps,p(ix,2)],varargin{:})-d(ix))/Deps; % Numerical
        dgrady = (feval(fd,[p(ix,1),p(ix,2)+Deps],varargin{:})-d(ix))/Deps; % gradient
        p(ix,:) = p(ix,:) - [d(ix).*dgradx,d(ix).*dgrady]; % Project back to boundary
    % 8. Termination criterion: All interior nodes move less than dptol (scaled)
        if max(sqrt(sum(deltat*Ftot(d<-Geps,:).^2,2))/h0) < dnptol
            break; 
        end
    end
end
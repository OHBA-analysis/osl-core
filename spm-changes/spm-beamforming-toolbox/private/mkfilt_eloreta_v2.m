function A=mkfilt_eloreta_v2(L,regu);
% makes spatial filter according to eLoreta 
% usage  A=mkfilt_eloreta_v2(L); or  A=mkfilt_eloreta_v2(L,regu);
%
% input L:  NxMxP leadfield tensor for N channels, M voxels, and 
%           P dipole directions. Typically P=3. (If you do MEG for 
%           a spherical volume conductor or reduce the rank, you must 
%           reduce L such that it has full rank for each voxel, such that,
%           e.g., P=2)
%       regu: optional regularization parameter (default is .05 corresponding 
%             to 5% of the average of the eigenvalues of some matrix to be inverted.) 
% 
% output A: NxMxP tensor of spatial filters. If x is the Nx1 data vector at time t. 
%           then A(:,m,p)'*x is the source activity at time t in voxel m in source direction
%           p. 
% 
% code implemented by Guido Nolte
% please cite
% 
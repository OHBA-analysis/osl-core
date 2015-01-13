ROInets.README.txt



%	Copyright 2015 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 263 $
%	$LastChangedDate: 2014-10-23 11:30:39 +0100 (Thu, 23 Oct 2014) $
%	Contact: giles.colclough 'at' eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 17-Mar-2014 23:39:39



ROInets - a Matlab package for performing leakage-robust network inference 
          between ROIs in MEG data. 

   The methodology used in this pipeline is set out in 
   Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., "A
   symmetric multivariate leakage correction for MEG connectomes,"
   NeuroImage [in revision]. 

This package provides tools for network analysis between regions of interest
in source-reconstructured MEG data, analysing amplitude correlations in 
band-limited power using a Gaussian Markov Random Field network model, after 
applying a symmetric multivariate orthogonalisation correction for source leakage. 



What you need for it to run:

FSL
FieldTrip
Matlab Stats + signal processing toolboxes - though you could probably change the code to make it work
QPAS mex files (optional)






What you need to get started:

Source-reconstructed resting-state MEG data (you can probably get it to work for task)
A set of ROIs or a spatial basis set, in the same space and 
  resolution as the MEG data, saved as a nifti






How to get started:

The +ROInets folder is a Matlab package. Don't change the name of the folder, 
and make sure the folder (not necessarily its contents) is visible on your path

View a brief summary of what each function does by typing `help ROInets.Contents'

The top-level function is ROInets.run_individual_network_analysis. 
View the helptext for this to view all the pipeline options. 

Look at the example file (`edit ROInets.example')

Most of the functions have an attempt at informative help text.







What paper to cite:
Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., "A
   symmetric multivariate leakage correction for MEG connectomes,"
   NeuroImage [in revision].








Where to ask for help, or report bugs:
get in touch with Giles at 
giles.colclough 'at' magd.ox.ac.uk
giles.colclough 'at' gmail.com






An overview of the pipeline:
- Select time region for analysis
- Band-pass filter the data
- Find a time-course for each ROI using PCA 
- Remove effects of leakage using an orthogonalisation process which finds the closest set of orthogonal vectors
- Find down-sampled power envelopes
- Perform network inference using partial correlation, L1 regularised using DP-glasso, optimised using cross-validation
- Convert correlations to z-statistics, by scaling relative to an empirical null, 
  generated to share the same temporal smoothness properties as the input data. 




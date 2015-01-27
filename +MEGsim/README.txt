osl_simulate_MEG_data.m code package
Released 16-Dec-2013
Included in osl 08-Apr-2014
$LastChangedBy: GilesColclough $
$Revision: 178 $
$LastChangedDate: 2013-12-17 11:37:58 +0000 (Tue, 17 Dec 2013) $

Giles Colclough
giles.colclough 'at' gmail.com
giles.colclough 'at' magd.ox.ac.uk


The main function is called osl_simulate_MEG_data.m. It is internally documented, and heavily commented.

Subfunctions are held in the package folder '+MEGsim'. Changing the name of this folder may cause runtime errors. Ensure the folder (but not necessarily its contents) is on the matlab path. 
FSL and OSL must also be set up correctly. See https://sites.google.com/site/ohbaosl/

A Matlab .mat file holds a covariance matrix describing noise measurements on an Electa Neuromag MEG system. This can be used to provide structured noise in the simulation. 

You can find a template MEEG object here:
https://sites.google.com/site/ohbaosl/practicals/practical-data



This software is released under the GNU General Public License v3 (or later), a copy of which is attached. 
It includes contributions by Mark Woolrich, Adam Baker and Tom Minka. 

Do please get in touch if you have questions or find bugs.

Giles

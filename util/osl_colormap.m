function map = osl_colormap(I)
  % OSL COLORMAPS
  % Perceptually uniform colormaps generated using COLORCET - see below
  %
  % VALID NAMES
  % - hot (FSL style red-yellow)
  % - cold (FSL style blue-light blue)
  % - electric (high-contrast blue)
  % - rwb (Classic red-white-blue)
  % - rkb (Red-black-blue)
  % - grey/gray (perceptually uniform greyscale)
  %
  % And there are also some primary colours that can be accessed
  % either by name or by index (passing in a number to this function instead of a string)
  %
  % - 1 : red
  % - 2 : blue
  % - 3 : green
  % - 4 : purple
  % - 5 : orange
  % - 6 : yellow
  % - 7 : brown
  % - 8 : pink
  %
  %
  %
  % Code is largely adapted from:
  %
  % Reference:
  % Peter Kovesi. Good Colour Maps: How to Design Them.
  % arXiv:1509.03700 [cs.GR] 2015
  % https://arxiv.org/abs/1509.03700

  % Copyright (c) 2013-2016 Peter Kovesi
  % Centre for Exploration Targeting
  % The University of Western Australia
  % peter.kovesi at uwa edu au
  % 
  % Permission is hereby granted, free of charge, to any person obtaining a copy
  % of this software and associated documentation files (the "Software"), to deal
  % in the Software without restriction, subject to the following conditions:
  % 
  % The above copyright notice and this permission notice shall be included in 
  % all copies or substantial portions of the Software.
  %
  % The Software is provided "as is", without warranty of any kind.

    N = 256;
    chromaK =   1;
    shift =   0;
    reverse =   0;
    diagnostics =   0;

    % Default parameters for colour map construction.  
    % Individual colour map specifications may override some of these.
    colourspace = 'LAB'; % Control points specified in CIELAB space
    sigma = 0;           % Default smoothing for colour map lightness equalisation.
    splineorder = 3;     % Default order of b-spline colour map curve.
    formula = 'CIE76';   % Default colour difference formula.
    W = [1 0 0];         % Colour difference weighting for Lightness,
                         % chroma and hue (default is to only use Lightness)
                         desc = '';
                         name = '';
                         attributeStr = '';
                         

   

    switch I   % Note all case expressions must be  case.

    %-----------------------------------------------------------------------------        
    %% Linear series

     case {'L1', 'grey', 'gray'}  % Grey 0 - 100 
      desc = 'Grey scale'; 
      attributeStr = 'linear';
      
      colpts = [  0 0 0      
                100 0 0];
      splineorder = 2;

    case {'L4a', 'hot'}
     desc = 'Red-Yellow heat colour map';
     attributeStr = 'linear';
     
     colourspace = 'RGB';
     splineorder = 2;       
     colpts = [
     1  0  0
     1  1  0];

  case {'L4b', 'cold'}
    desc = 'Blue to light blue';
    attributeStr = 'linear';
    
    colourspace = 'RGB';
    splineorder = 2;       
    colpts = [
    0  0.15  1
    0  1  1];

  case {'L6','electric'}
    desc = 'Blue shades running vertically up the blue edge of CIELAB space';
    attributeStr = 'linear';
    
    colpts = [ 5  31 -45
    15 50 -66
    25 65 -90
    35 70 -100
    45 45 -85
    55  20 -70
    65  0 -53
    75 -22 -37
    85 -38 -20
    95 -25 -3]; 

  case {1,'red'}                  
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([228,26,28]);
  case {2,'blue'}                  
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([55,126,184]);
  case {3,'green'}                 
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([77,175,74]);
  case {4,'purple'}                 
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([152,78,163]);
  case {5,'orange'}
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([255,127,0]);
  case {6,'yellow'}
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([255,255,51]);
  case {7,'brown'}
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([166,86,40]);
  case {8,'pink'}
      [attributeStr,splineorder,colourspace,colpts] = map_sequence([247,129,191]);

      case {'L2', 'REDUCEDGREY'} %  Grey 10 - 95
       desc = ['Grey scale with slightly reduced contrast to '...
       'avoid display saturation problems'];
       attributeStr = 'linear';
       
       colpts = [10 0 0
       95 0 0];
       splineorder = 2;
       


       case 'L5'
         desc = 'Colour Map along the green edge of CIELAB space';
         attributeStr = 'linear';
         
         colpts = [ 5 -9  5
         15 -23 20
         25 -31 31
         35 -39 39
         45 -47 47
         55 -55 55
         65 -63 63
         75 -71 71
         85 -79 79
         95 -38 90]; 
         
         % Heat map from straight line segments from black to red to yellow to white
       case {'L3A'}
         desc = 'Black-Red-Yellow-White heat colour map';
         attributeStr = 'linear';
         
         colourspace = 'RGB';
         splineorder = 2;       
         colpts = [0 0 0
         1 0 0
         1 1 0
         1 1 1 ];



           % Colour map as above but with spline order set to 3 and additional
           % control points to get a similar colour map path but with the corners
           % at red and yellow rounded off.  Works well
         case {'L3', 'HEAT', 'FIRE'}
           desc = 'Black-Red-Yellow-White heat colour map';
           attributeStr = 'linear';
           
           colourspace = 'RGB';
           splineorder = 3;       
           colpts = [0 0 0
           .85 0 0
           1 .15 0
           1 .85 0
           1 1 .15
           1 1 1 ];      

         case {'L4', 'HEATYELLOW'}
           desc = 'Black-Red-Yellow heat colour map';
           attributeStr = 'linear';
           
           colourspace = 'RGB';
           splineorder = 3;       
           colpts = [0 0 0
           .85 0 0
           1 .15 0
           1  1  0];

  case 'L7'
    desc = 'Blue-Pink-Light Pink colour map';
    attributeStr = 'linear';
    
    colpts = [ 5 29 -43
    15 48 -66
    25 64 -89
    35 73 -100
    45 81 -86
    55 90 -69
    65 83 -53
    75 56 -36
    85 32 -22
    95 10 -7];

  case {'L8', 'BMY'}
    desc = 'Blue-Magenta-Yellow highly saturated colour map';
    attributeStr = 'linear';
    
    colpts = [10 ch2ab( 55,-58)
    20 ch2ab( 75,-58)
    30 ch2ab( 75,-40)
    40 ch2ab( 73,-20)
    50 ch2ab( 75,  0)                 
    60 ch2ab( 70, 30)
    70 ch2ab( 65, 60)
    80 ch2ab( 75, 80)
    95 ch2ab( 80, 105)];


     case {'L9', 'BGYW'}   % Blue to yellow section of R1 with short
                           % extensions at each end 
                           desc = 'Blue-green-yellow colour map';
                           attributeStr = 'linear';
                           
                           colpts = [15  49 -64
                           35  67 -100
                           45 -12 -29
                           60 -55 60
                           85 -10 80
                           95 -17 50
                           100 0  0];

                         case {'L10', 'GEOGRAPHIC'}
                          desc = ['A "geographical" colour map.  '...
                          'Best used with relief shading'];
                          attributeStr = 'linear';
                          
      colpts = [60 ch2ab(20, 180)   % pale blue green
      65 ch2ab(30, 135)
      70 ch2ab(35, 75)
      75 ch2ab(45, 85)
      80 ch2ab(22, 90)        
      85 0 0   ];
      
     case {'L11', 'GEOGRAPHIC2'}  %  Lighter version of L10 with a bit more chroma 
      desc = ['A lighter "geographical" colour map.  '...
      'Best used with relief shading'];
      attributeStr = 'linear';
      
      colpts = [65 ch2ab(50, 135)   % pale blue green
      75 ch2ab(45, 75)
      80 ch2ab(45, 85)
      85 ch2ab(22, 90)        
      90 0 0   ];

    case {'L12', 'DEPTH'}
      desc =  'A "water depth" colour map';
      attributeStr = 'linear';
      
      colpts = [95 0 0
      80 ch2ab(20, -95)
      70 ch2ab(25, -95)
      60 ch2ab(25, -95)
      50 ch2ab(35, -95)];

      
     % The following three colour maps are for ternary images, eg Landsat images
     % and radiometric images.  These colours form the nominal red, green and
     % blue 'basis colours' that are used to form the composite image.  They are
     % designed so that they, and their secondary colours, have nearly the same
     % lightness levels and comparable chroma.  This provides consistent feature
     % salience no matter what channel-colour assignment is made.  The colour
     % maps are specified as straight lines in RGB space.  For their derivation
     % see
     % http://peterkovesi.com/projects/colourmaps/ColourmapTheory/index.html#ternary
     
   case {'L13', 'REDTERNARY'}
    desc = 'red colour map for ternary images';
    attributeStr = 'linear';
    
    colourspace = 'RGB';
    colpts = [0.00 0.00 0.00
    0.90 0.17 0.00];

    splineorder = 2; 

  case {'L14', 'GREENTERNARY'}
    desc = 'green colour map for ternary images';
    attributeStr = 'linear';
    
    colourspace = 'RGB';
    colpts = [0.00 0.00 0.00
    0.00 0.50 0.00];

    splineorder = 2; 

  case {'L15', 'BLUETERNARY'}
    desc = 'blue colour map for ternary images';
    attributeStr = 'linear';
    
    colourspace = 'RGB';
    colpts = [0.00 0.00 0.00
    0.10 0.33 1.00];

    splineorder = 2; 







     %--------------------------------------------------------------------------------
     %% Diverging colour maps
     
     % Note that on these colour maps often we do not go to full white but use a
     % lightness value of 95. This helps avoid saturation problems on monitors.
     % A lightness smoothing sigma of 5 to 7 is used to avoid generating a false
     % feature at the white point in the middle.  Note however, this does create
     % a small perceptual contrast blind spot at the middle.

   case {'D1', 'rwb'}
    desc = ['Classic diverging blue - white - red colour map. End colours are' ...
    ' matched in lightness and chroma'];
    attributeStr = 'diverging';
    
    colpts = [40  ch2ab(83,-64)
    95  0   0
    40  ch2ab(83, 39)];   
    sigma = 7;
    splineorder = 2; 


     case 'D1A'   % Attempt to improve metric property of D1
      desc = 'Diverging blue - white - red colour map';
      attributeStr = 'diverging';
      
      colpts = [40  ch2ab(83,-64)
      65  ch2ab(70, -28)
      95  0   0
      95  0   0
%                65  ch2ab(60, 13)
65  ch2ab(70, 75)
40  ch2ab(83, 39)];   
sigma = 7;
splineorder = 3; 


case 'D2'
  desc = 'Diverging green - white - violet colour map';
  attributeStr = 'diverging';
  
  colpts = [55 -50  55
  95   0   0
  55  60 -55];  
  sigma = 7;
  splineorder = 2; 

case 'D3'
  desc = 'Diverging green - white - red colour map';
  attributeStr = 'diverging';
  
  colpts = [55 -50 55
  95   0  0
  55  63 39];  
  sigma = 7;
  splineorder = 2;  

case {'D1', 'rkb'}
  desc = 'Diverging blue - black - red colour map';
  attributeStr = 'diverging';
  
  colpts = [55 ch2ab(70, -76)
  10  0   0
  55 ch2ab(70, 35)];   
  sigma = 7;
  splineorder = 2; 

case 'D5'
  desc = 'Diverging green - black - red colour map';
  attributeStr = 'diverging';
  
  colpts = [60 ch2ab(80, 134)
  10  0   0
  60 ch2ab(80, 40)];   
  sigma = 7;
  splineorder = 2; 

case 'D6'
  desc = 'Diverging blue - black - yellow colour map';
  attributeStr = 'diverging';
  
  colpts = [60 ch2ab(60, -85)
  10  0   0
  60 ch2ab(60, 85)];   
  sigma = 7;
  splineorder = 2;        

     case {'D7', 'DIVBGY'}  % Linear diverging  blue - grey - yellow.  Works well
      desc = ['Linear-diverging blue - grey - yellow colour map. This kind' ...
      ' of diverging map has no perceptual dead spot at the centre'];
      attributeStr = 'diverging-linear';
      
      colpts = [30 ch2ab(89, -59)
      60 0 0                 
      90 ch2ab(89,96)];
      splineorder = 2;
      
     case 'D7B'  % Linear diverging blue - grey - yellow.
                 % Similar to 'D7' but with slight curves in the path to
                 % slightly increase the chroma at the 1/4 and 3/4 locations
                 % in the map.
                 desc = ['Linear-diverging blue - grey - yellow colour map. This kind' ...
                 ' of diverging map has no perceptual dead spot at the centre'];                 
                 attributeStr = 'diverging-linear';
                 

                 h1 = -59; h2 = 95;
                 mh = (h1+h2)/2;
                 dh = 10;

                 colpts = [30 ch2ab(88, h1)
                 48 ch2ab(40, h1+dh)
                 60 0 0       
                 60 0 0       
                 72 ch2ab(40, h2-dh)
                 90 ch2ab(88, h2)];
                 splineorder = 3;

     case 'D8' % Linear diverging  blue - grey - red
      desc = 'Linear-diverging blue - grey - red';
      attributeStr = 'diverging-linear';
      
      colpts = [30 ch2ab(105, -58)
      42.5 0 0                 
      55 ch2ab(105,41)];
      splineorder = 2;
      
     case 'D9'  % Lightened version of D1 for relief shading - Good.
      desc = ['Diverging low contrast blue - white - red colour map.  Good in' ...
      ' conjunction with relief shading'];
      attributeStr = 'diverging';
      
      colpts = [55  ch2ab(73,-74)
      98  0   0
      55  ch2ab(73, 39)];   
      sigma = 2;   % Less smoothing needed for low contrast
      splineorder = 2; 
      
     case 'D10' % low contrast diverging map for when you want to use
                % relief shading 
                desc = ['Diverging low contrast cyan - white - magenta colour map.  Good in' ...
                ' conjunction with relief shading'];
                attributeStr = 'diverging';
                
                colpts = [80 ch2ab(44, -135)
                100  0   0
                80 ch2ab(44, -30)];   
      sigma = 0;   % No smoothing needed for lightness range of 20
      splineorder = 2; 
      W = [1 1 1];
      
     case 'D11' % Constant lightness diverging map for when you want to use
                % relief shading  ? Perhaps lighten the grey so it is not quite
                % isoluminant ?
                desc = 'Diverging isoluminant lightblue - lightgrey - orange colour map';
                attributeStr = 'diverging-isoluminant';
                
                colpts = [70 ch2ab(50, -105)
                70  0   0
                70 ch2ab(50, 45)];   
                sigma = 7;
                splineorder = 2; 
                W = [1 1 1];

     case 'D12' % Constant lightness diverging map for when you want to use
                % relief shading  ? Perhaps lighten the grey so it is not quite
                % isoluminant ?
                desc = 'Diverging isoluminant lightblue - lightgrey - pink colour map';
                attributeStr = 'diverging-isoluminant';
                
                colpts = [75 ch2ab(46, -122)
                75  0   0
                75 ch2ab(46, -30)];   
                sigma = 7;
                splineorder = 2; 
                W = [1 1 1];


     %-------------------------------------------------------------------------      
     %% Cyclic colour maps

     case 'C1' % I think this is my best zigzag style cyclic map - Good!
               % Control points are placed so that lightness steps up and
               % down are equalised.  Additional intermediate points are
               % placed to try to even up the 'spread' of the key colours.
               desc = ['Cyclic: magenta - red - yellow - blue. Alows four' ...
               ' orientations/phase angles to be visulaised.'];
               attributeStr = 'cyclic';
               
               mag = [75 60 -40];
               yel = [75 0 77];
               blu = [35  70 -100];
               red = [35 60 48];
               colpts = [mag
               55 70 0 
               red
               55 35 62
               yel
               50 -20 -30
               blu
               55 45 -67
               mag];

               sigma = 7;
       splineorder = 2;  % linear path 

     case {'C2', 'PHASE4'}  % A big diamond across the gamut.  Really good!  Incorporates two
               % extra cotnrol points around blue to extend the width of that
               % segment slightly.
               desc = ['Cyclic: magenta - yellow - green - blue. Alows four' ...
               ' orientations/phase angles to be visulaised.  Really good!'];
               attributeStr = 'cyclic';
               
               colpts = [62.5  83 -54
               80 20 25
               95 -20 90
               62.5 -65 62   
               42 10 -50
               30 75 -103                
               48 70 -80  
               62.5  83 -54];
               sigma = 7;
               splineorder = 2;

     case 'C3'   % red-white-blue-black-red allows quadrants to be identified
      desc = ['Cyclic: red - white - blue - black - red.  Allows four' ...
      ' orientations/phase angles to be visulaised.'];
      attributeStr = 'cyclic';
      
      colpts = [50 ch2ab(85, 39)
      85 0 0
      50 ch2ab(85, -70)
      15 0 0
      50 ch2ab(85, 39)];
      sigma = 7;
      splineorder = 2;
      
     case {'C4', 'PHASE2'}   % white-red-white-blue-white Works nicely
      desc = ['Cyclic: white - red - white - blue. Good if you just want to' ...
      ' visualise +ve and -ve phase'];  
      attributeStr = 'cyclic';
      
      colpts = [90 0 0
      40 65 56
      90 0 0
      40 31 -80
      90 0 0];
      sigma = 7;
      splineorder = 2;      

     case {'C5', 'CYCLICGREY'}   % Cyclic greyscale  Works well
      desc = 'Cyclic: greyscale'; 
      attributeStr = 'cyclic';
      
      colpts = [50 0 0
      85 0 0
      15 0 0
      50 0 0];
      sigma = 7;
      splineorder = 2;
      
     case 'C6'  % Circle at 67  - sort of ok but a bit flouro
      desc = 'Cyclic: isoluminant';
      attributeStr = 'cyclic-isoluminant colour circle';
      
      chr = 42;
      ang = 124;
      colpts = [67  ch2ab(chr,  ang-90)
      67  ch2ab(chr,  ang)
      67  ch2ab(chr,  ang+90)
      67  ch2ab(chr,  ang+180)
      67  ch2ab(chr,  ang-90)];
      W = [1 1 1];
      
     case 'C7'  % Elliptical path - ok
      desc = 'Cyclic: low contrast colour ellipse';
      attributeStr = 'cyclic';
      
      ang = 112;
      colpts = [70    ch2ab(46,  ang-90)
      90    ch2ab(82,  ang)
      70    ch2ab(46,  ang+90)
      50    ch2ab(82,  ang+180)
      70    ch2ab(46,  ang-90)];
      W = [1 1 1];
      
     case 'C8' % Elliptical path.  Greater range of lightness values and
               % slightly more saturated colours.  Seems to work however I do
               % not find the colour sequence that attractive. This is a
               % constraint of the gamut.
               desc = 'Cyclic: colour ellipse';
               attributeStr = 'cyclic';
               
               ang = 124;
               colpts = [60    ch2ab(40,  ang-90)
               100   ch2ab(98,  ang)
               60    ch2ab(40,  ang+90)
               20    ch2ab(98,  ang+180)
               60    ch2ab(40,  ang-90)];
               W = [1 1 1];
               sigma = 7;                        

     case 'C9' % Variation of C1. Perceptually this is good. Excellent balance
               % of colours in the quadrants but the colour mix is not to my
               % taste.  Don't like the green.  The red-green transition
               % clashes
               desc = ['Cyclic: blue - green - red - magenta. Allows four' ...
               ' orientations/phase angles to be visulaised.']; 
               attributeStr = 'cyclic';
               
               blu = [35  70 -100];
               colpts = [blu
               70 -70 64
               35 65 50
               70 75 -46
               blu        ];      
               sigma = 7;
       splineorder = 2;  % linear path 


     %-----------------------------------------------------------------------------    
     %%  Rainbow style colour maps

     case {'R1', 'RAINBOW'}   % Reasonable rainbow colour map after it has been
                              % fixed by equalisecolourmap.
                              desc = ['The least worst rainbow colour map I can devise.  Note there are' ...
                              ' small perceptual blind spots at yellow and red'];
                              attributeStr = 'rainbow';
                              
                              colpts = [35 63 -98
                              45 -14 -30
                              60 -55 60
                              85 0 80
                              55 60 62
                              75 55 -35];      
                              sigma = 7;
       splineorder = 2;  % linear path 

     case {'R2', 'RAINBOW2'}   % Similar to R1 but with the colour map finishing
                               % at red rather than continuing onto pink.
                               desc = ['Reasonable rainbow colour map from blue to red.  Note there is' ...
                               ' a small perceptual blind spot at yellow'];
                               attributeStr = 'rainbow';
                               
                               colpts = [35 60 -98
                               45 -15 -30
                               60 -55 60
                               85 0 78
                               55 73 68];
                               sigma = 5;
       splineorder = 2;  % linear path 
       
       
     case {'R3', 'RAINBOW3'}   % Diverging rainbow.  The blue and red points are matched in
                               % lightness and chroma as are the green and magenta points 
                               desc = ['Diverging-rainbow colourmap. Yellow is the central reference' ...
                               ' colour. The blue and red end points are matched in lightness' ...
                               ' and chroma'];
                               attributeStr = 'diverging-rainbow';
                               
                               colpts = [45 39 -83
                               52 -23 -23
                               60 -55 55
                               85 -2 85
                               60 74 -17
                               45 70 59];

                               sigma = 5;
       splineorder = 2;  % linear path 
       
       

     %-----------------------------------------------------------------------------    
     %%  Isoluminant colour maps       

   case 'I1'
    desc = ['Isoluminant blue to green to orange at lightness 70.  '...
    'Poor on its own but works well with relief shading'];
    attributeStr = 'isoluminant';
    
    colpts = [70 ch2ab(40, -115)
    70 ch2ab(50, 160)
    70 ch2ab(50,  90)
    70 ch2ab(50,  45)];
    W = [1 1 1];

     case 'I2'  % Adaptation of I1 shifted to 80 from 70
      desc = ['Isoluminant blue to green to orange at lightness 80.  '...
      'Poor on its own but works well with relief shading'];
      attributeStr = 'isoluminant';
      
      colpts = [80 ch2ab(36, -115)
      80 ch2ab(50, 160)
      80 ch2ab(50,  90)
      80 ch2ab(46,  55)];
      W = [1 1 1];
      
    case 'I3'
      desc = ['Isoluminant blue to pink at lightness 70.  '...
      'Poor on its own but works well with relief shading'];
      attributeStr = 'isoluminant';
      
      colpts = [70 ch2ab(40, -125)
      70 ch2ab(40, -80)
      70 ch2ab(40, -40)
      70 ch2ab(50,  0)];
      W = [1 1 1];      
      

      
     %-----------------------------------------------------------------------------    
     %%  Experimental colour maps and colour maps that illustrate some design principles


   case 'X1'
    desc = ['Two linear segments with different slope to illustrate importance' ...
    ' of lightness gradient.'];
    attributeStr = 'linear-lightnessnormalised';
    
    colpts = [30 ch2ab(102, -54)
    40  0   0
    90  ch2ab(90, 95)];
    W = [1 0 0];
    splineorder = 2;

  case 'X2'
    desc = ['Two linear segments with different slope to illustrate importance' ...
    ' of lightness gradient.'];
    attributeStr = 'linear-CIE76normalised';
    
    colpts = [30 ch2ab(102, -54)
    40  0   0
    90  ch2ab(90, 95)];
    W = [1 1 1];
    splineorder = 2;      


     case 'X3' % Constant lightness 'v' path to test unimportance of having a smooth
               % path in hue/chroma.  Slight 'feature' at the red corner (Seems more
               % important on poor monitors)
               attributeStr = 'isoluminant-HueChromaSlopeDiscontinuity';
               
               colpts = [50 17 -78
               50 77 57
               50 -48 50];
      splineorder = 2;  % linear path      
      W = [1 1 1];      
      


     % A set of isoluminant colour maps only varying in saturation to test
     % the importance of saturation (not much) Colour Maps are linear with
     % a reversal to test importance of continuity.  

     case 'X10'  % Isoluminant 50 only varying in saturation
      attributeStr = 'isoluminant';
      
      colpts = [50 0 0
      50 77 64
      50 0 0];
      splineorder = 2;
      sigma = 0;
      W = [1 1 1];
      
     case 'X11'  % Isoluminant 50 only varying in saturation
      attributeStr = 'isoluminant';
      
      colpts = [50 0 0
      50 0 -56
      50 0 0];
      splineorder = 2;
      sigma = 0;
      W = [1 1 1];      
      
     case 'X12'  % Isoluminant 90 only varying in saturation
      attributeStr = 'isoluminant';
      
      colpts = [90 0 0
      90 -76 80
      90 0 0];
      splineorder = 2;
      sigma = 0;
      W = [1 1 1];      

      
    % Difference in CIE76 and CIEDE2000 in chroma      

    case 'X13'  % Isoluminant 55 only varying in chroma. CIEDE76
      attributeStr= 'isoluminant-CIE76';
      
      colpts = [55 0 0
      55 80 67];
      splineorder = 2;
      W = [1 1 1];      
      formula = 'CIE76';
      
     case 'X14'  % Same as X13 but using CIEDE2000
      attributeStr= 'isoluminant-CIEDE2000';
      
      colpts = [55 0 0
      55 80 67];
      splineorder = 2;
      W = [1 1 1];      
      formula = 'CIEDE2000';      
      
     case 'X15'  % Grey 0 - 100. Same as No 1 but with CIEDE2000
      desc = 'Grey scale'; 
      attributeStr= 'linear-CIEDE2000';
      
      colpts = [  0 0 0      
      100 0 0];
      splineorder = 2;
      formula = 'CIEDE2000';   

     case 'X16'  % Isoluminant 30 only varying in chroma
      attributeStr= 'isoluminant';
      
      colpts = [30 0 0
      30 77 -106];
      splineorder = 2;
      W = [1 1 1];      

      
     case 'X21' % Blue to yellow section of rainbow map R1 for illustrating
                % colour ordering issues
                attributeStr= 'rainbow-section1';
                
                colpts = [35 60 -100
                45 -15 -30
                60 -55 60
                85 0 80];
       splineorder = 2;  % linear path 

     case 'X22' % Red to yellow section of rainbow map R1 for illustrating
                % colour ordering issues
                attributeStr= 'rainbow-section2';
                
                colpts = [55 70 65
                85 0 80];
       splineorder = 2;  % linear path 

     case 'X23' % Red to pink section of rainbow map R1 for illustrating
                % colour ordering issues
                attributeStr= 'rainbow-section3';                
                
                colpts = [55 70 65
                75 55 -35];      
       splineorder = 2;  % linear path 

       
     case 'XD1'  % Same as D1 but with no smoothing
      desc = 'Diverging blue-white-red colour map';
      attributeStr = 'diverging';
      
      colpts = [40  ch2ab(83,-64)
      95  0   0
      40  ch2ab(83, 39)];   
      sigma = 0;
      splineorder = 2; 
      

    case 'X30'
      desc = 'red - green - blue interpolated in rgb';
      attributeStr = 'linear';
      
      colourspace = 'RGB';
      colpts = [1.00 0.00 0.00
      0.00 1.00 0.00
      0.00 0.00 1.00];
      sigma = 0;
      W = [0 0 0];
      splineorder = 2; 
      
    case 'X31'
      desc = 'red - green - blue interpolated in CIELAB';
      attributeStr = 'linear';
      
      colourspace = 'LAB';
      colpts = [53  80   67
      88 -86   83
      32  79 -108];
      sigma = 0;
      W = [0 0 0];
      splineorder = 2;       

     case 'XD7A'  % Linear diverging blue - magenta- grey - orange - yellow.
                 % Modified from 'D7' to have a double arch shaped path in an attempt
                 % to improve its Metric properties.  Also starts at lightness
                 % of 40 rather than 30.  The centre grey region is a bit too
                 % prominant and overall the map is perhaps a bit too 'bright'
                 attributeStr = 'diverging-linear';
                 
                 colpts = [40 ch2ab(88, -64)
                 55 ch2ab(70, -30)
                 64 ch2ab(2.5, -72.5)
                 65 0 0       
                 66 ch2ab(2.5, 107.5)
                 75 ch2ab(70, 70)
                 90 ch2ab(88,100)];
                 splineorder = 3;


     case 'XC15B'   % Francesca Samsel's beautiful c15b blue-white-green
                    % diverging map.  Map is reproduced in its original form
                    % with no equalisation of lightness gradient.  (It is
                    % close to being equal as originally designed.)
                    desc = 'Francesca Samsel''s c15b blue-white-green diverging map';
                    attributeStr = 'diverging';
                    
                    colourspace = 'RGB';
                    colpts = [0.231373 0.247059 0.329412
                    0.266667 0.305882 0.411765
                    0.286275 0.368627 0.478431
                    0.301961 0.439216 0.549020
                    0.309804 0.521569 0.619608
                    0.380392 0.631373 0.690196
                    0.454902 0.745098 0.760784
                    0.541176 0.831373 0.803922
                    0.631373 0.901961 0.843137
                    0.768627 0.960784 0.894118
                    0.901961 1.000000 0.949020
                    0.768627 0.960784 0.835294
                    0.635294 0.909804 0.698039
                    0.552941 0.850980 0.576471
                    0.490196 0.780392 0.466667
                    0.447059 0.701961 0.384314
                    0.407843 0.611765 0.305882
                    0.360784 0.509804 0.231373
                    0.305882 0.400000 0.160784
                    0.231373 0.278431 0.098039
                    0.141176 0.149020 0.043137];

                    sigma = 0;  
      W = [0 0 0];  % Set to zero to reproduce original map exactly
      splineorder = 2;             

     case 'XC15BM'   % Francesca's map modified slightly with additional
                     % starting control point to create symmetric lightness
                     % profile. Alternatively, remove the near black colour at
                     % the end.  Map is equalised for lightness gradient.  A
                     % really nice map!
                     desc = 'Francesca Samsel''s c15b blue-white-green diverging map';
                     attributeStr = 'diverging';
                     
                     colourspace = 'RGB';
      colpts = [%0.120752 0.138784 0.214192 % additional point to match lightness at end
      0.231373 0.247059 0.329412
      0.266667 0.305882 0.411765
      0.286275 0.368627 0.478431
      0.301961 0.439216 0.549020
      0.309804 0.521569 0.619608
      0.380392 0.631373 0.690196
      0.454902 0.745098 0.760784
      0.541176 0.831373 0.803922
      0.631373 0.901961 0.843137
      0.768627 0.960784 0.894118
      0.901961 1.000000 0.949020
      0.768627 0.960784 0.835294
      0.635294 0.909804 0.698039
      0.552941 0.850980 0.576471
      0.490196 0.780392 0.466667
      0.447059 0.701961 0.384314
      0.407843 0.611765 0.305882
      0.360784 0.509804 0.231373
      0.305882 0.400000 0.160784
%                0.231373 0.278431 0.098039
                  0.2218    0.2690    0.0890 % tweeked to match start lightness
               % 0.141176 0.149020 0.043137
               ];

               sigma = 7;  
               W = [1 0 0];
               splineorder = 2;             


     case 'XD7C'  % Linear diverging  green - grey - yellow
                  % Reasonable, perhaps easier on the eye than D7

                  attributeStr = 'diverging-linear';
                  
                  rad = 65;
                  colpts = [40 ch2ab(rad, 136)
                  65 0 0                 
                  90 ch2ab(rad, 95)];
                  splineorder = 2;


                end


    % Adjust chroma/saturation but only if colourspace is LAB
    if strcmpi(colourspace, 'LAB')
      colpts(:,2:3) = chromaK * colpts(:,2:3);
    end
    
    % Colour map path is formed via a b-spline in the specified colour space
    Npts = size(colpts,1);
    
    if Npts < 2
      error('Number of input points must be 2 or more')
    elseif Npts < splineorder
      splineorder = Npts;
      fprintf('Warning: Spline order is greater than number of data points\n')
      fprintf('Reducing order of spline to %d\n', Npts)
    end
    
    % Rely on the attribute string to identify if colour map is cyclic.  We may
    % want to construct a colour map that has identical endpoints but do not
    % necessarily want continuity in the slope of the colour map path.
    if strfind(attributeStr, 'cyclic')
      cyclic = 1;
      labspline = pbspline(colpts', splineorder, N);
    else
      cyclic = 0;
      labspline = bbspline(colpts', splineorder, N);
    end    
    
    % Apply contrast equalisation with required parameters. Note that sigma is
    % normalised with respect to a colour map of length 256 so that if a short
    % colour map is specified the smoothing that is applied is adjusted to suit.
    sigma = sigma*N/256;
    map = equalisecolourmap(colourspace, labspline', formula,...
      W, sigma, cyclic, diagnostics);    

    % If specified apply a cyclic shift to the colour map
    if shift
      if isempty(strfind(attributeStr, 'cyclic'))
        fprintf('Warning: Colour map shifting being applied to a non-cyclic map\n');
      end
      map = circshift(map, round(N*shift));    
    end
    
    if reverse
     map = flipud(map);
   end    

    % Compute mean chroma of colour map
    lab = rgb2lab(map);
    meanchroma = sum(sqrt(sum(lab(:,2:3).^2, 2)))/N;
    
    % Construct lightness range description
    if strcmpi(colourspace, 'LAB')  % Use the control points
      L = colpts(:,1);
    else  % For RGB use the converted CIELAB values
      L = round(lab(:,1));
    end
    minL = min(L);
    maxL = max(L);
    
    if minL == maxL     % Isoluminant
      LStr = sprintf('%d', minL);

    elseif L(1) == maxL && ...
     (~isempty(strfind(attributeStr, 'diverging')) ||...
      ~isempty(strfind(attributeStr, 'linear')))
     LStr = sprintf('%d-%d', maxL, minL);

   else          
    LStr = sprintf('%d-%d', minL, maxL);
  end

%------------------------------------------------------------------
% Conversion from (chroma, hue angle) description to (a*, b*) coords

function ab = ch2ab(chroma, angle_degrees)

  theta = angle_degrees/180*pi;
  ab = chroma*[cos(theta) sin(theta)];



% BBSPLINE - Basic B-spline
%
% Usage:  S = bbspline(P, k, N)
% 
% Arguments:   P - [dim x Npts] array of control points
%              k - order of spline (>= 2). 
%                  k = 2: Linear
%                  k = 3: Quadratic, etc
%              N - Optional number of points to evaluate along
%                  spline. Defaults to 100.
%
% Returns:     S - spline curve  [dim x N] spline points
%
% See also: PBSPLINE

% PK Jan 2014
%    Nov 2015  Made basis calculation slightly less wasteful

function S = bbspline(P, k, N)

  if ~exist('N', 'var'), N = 100; end

    [dim, np1] = size(P);
    n = np1-1;

    assert(k >= 2, 'Spline order must be 2 or greater');    
    assert(np1 >= k, 'No of control points must be >= k');
    assert(N >= 2, 'Spline must be evaluated at 2 or more points');

    % Set up open uniform knot vector from 0 - 1.  
    % There are k repeated knots at each end.
    ti = 0:(k+n - 2*(k-1));    
    ti = ti/ti(end);
    ti = [repmat(ti(1), 1, k-1), ti, repmat(ti(end), 1, k-1)];

    nK = length(ti);
    
    % Generate values of t that the spline will be evaluated at
    dt = (ti(end)-ti(1))/(N-1);
    t = ti(1):dt:ti(end);
    
    % Build complete array of basis functions.  We maintain two
    % arrays, one storing the basis functions at the current level of
    % recursion, and one storing the basis functions from the previous
    % level of recursion
    B = cell(1,nK-1);
    Blast = cell(1,nK-1);
    
    % 1st level of recursive construction
    for i = 1:nK-1
     Blast{i} = t >= ti(i) & t < ti(i+1)  & ti(i) < ti(i+1); 
   end

    % Subsequent levels of recursive basis construction.  Note the logic to
    % handle repeated knot values where ti(i) == ti(i+1)
    for ki = 2:k    
      for i = 1:nK-ki

        if (ti(i+ki-1) - ti(i)) < eps
          V1 = 0;
        else
          V1 = (t - ti(i))/(ti(i+ki-1) - ti(i)) .* Blast{i};
        end

        if (ti(i+ki) - ti(i+1)) < eps
          V2 = 0;
        else
          V2 = (ti(i+ki) - t)/(ti(i+ki) - ti(i+1)) .* Blast{i+1};
        end

        B{i} = V1 + V2;

%           This is the ideal equation that the code above implements            
%            B{i,ki} = (t - ti(i))/(ti(i+ki-1) - ti(i)) .* B{i,ki-1} + ...
%                      (ti(i+ki) - t)/(ti(i+ki) - ti(i+1)) .* B{i+1,ki-1};

end

        % Swap B and Blast, but only if this is not the last iteration
        if ki < k
          tmp = Blast;
          Blast = B;
          B = tmp;
        end
      end

    % Apply basis functions to the control points
    S = zeros(dim, length(t));
    for d = 1:dim
      for i = 1:np1
        S(d,:) = S(d,:) + P(d,i)*B{i};
      end
    end
    
    % Set the last point of the spline. This is not evaluated by the code above
    % because the basis functions are defined from ti(i) <= t < ti(i+1)
    S(:,end) = P(:,end);

    function newrgbmap = equalisecolourmap(rgblab, map, formula, W, sigma, cyclic, diagnostics)

      if ~exist('formula', 'var'), formula = 'CIE76'; end
        if ~exist('W', 'var'), W = [1 0 0]; end
          if ~exist('sigma', 'var'), sigma = 0; end
            if ~exist('cyclic', 'var'), cyclic = 0; end
              if ~exist('diagnostics', 'var'), diagnostics = 0; end

    N = size(map,1);  % No of colour map entries

    if N/sigma < 25
      warning(['It is not recommended that sigma be larger than 1/25 of' ...
       ' colour map length'])
    end

    if strcmpi(rgblab, 'RGB') && (max(map(:)) > 1.01  || min(map(:) < -0.01))
      error('If map is RGB values should be in the range 0-1')
    elseif strcmpi(rgblab, 'LAB') && max(abs(map(:))) < 10
      error('If map is LAB magnitude of values are expected to be > 10')
    end
    
    % Use D65 whitepoint to match typical monitors.
    wp = whitepoint('D65');
    
    % If input is RGB convert colour map to Lab. 
    if strcmpi(rgblab, 'RGB')
      rgbmap = map;
      labmap = applycform(map, makecform('srgb2lab', 'AdaptedWhitePoint', wp));
      L = labmap(:,1);
      a = labmap(:,2);
      b = labmap(:,3);
    elseif strcmpi(rgblab, 'LAB')
      labmap = map;
      rgbmap = applycform(map, makecform('lab2srgb', 'AdaptedWhitePoint', wp));
      L = map(:,1);
      a = map(:,2);
      b = map(:,3);
    else
      error('Input must be RGB or LAB')
    end

    % The following section of code computes the locations to interpolate into
    % the colour map in order to achieve equal steps of perceptual contrast.
    % The process is repeated recursively on its own output. This helps overcome
    % the approximations induced by using linear interpolation to estimate the
    % locations of equal perceptual contrast. This is mainly an issue for
    % colour maps with only a few entries.

    for iter = 1:3
        % Compute perceptual colour difference values along the colour map using
        % the chosen formula and weighting vector.
        if strcmpi(formula, 'CIE76')
          deltaE = cie76(L, a, b, W);
        elseif strcmpi(formula, 'CIEDE2000')
          if ~osl_util.isfile('deltaE2000')
            fprintf('You need Gaurav Sharma''s function deltaE2000.m\n');
            fprintf('available at:\n');
            fprintf('http://www.ece.rochester.edu/~gsharma/ciede2000/\n');
            newrgbmap = [];
            return
          end

          deltaE = ciede2000(L, a, b, W);
        else
          error('Unknown colour difference formula')
        end
        
        % Form cumulative sum of of delta E values.  However, first ensure all
        % values are larger than 0.001 to ensure the cumulative sum always
        % increases.
        deltaE(deltaE < 0.001) = 0.001;
        cumdE = cumsum(deltaE);
        
        % Form an array of equal steps in cumulative contrast change.
        equicumdE =  (0:(N-1))'/(N-1) * (cumdE(end)-cumdE(1)) + cumdE(1);
        
        % Solve for the locations that would give equal Delta E values.
        method = 'linear';
        newN = interp1(cumdE, (1:N)', equicumdE, method, 'extrap');
        
        % newN now represents the locations where we want to interpolate into the
        % colour map to obtain constant perceptual contrast
        L = interp1((1:N)', L, newN, method, 'extrap');  
        a = interp1((1:N)', a, newN, method, 'extrap');
        b = interp1((1:N)', b, newN, method, 'extrap');

        % Record initial colour differences for evaluation at the end
        if iter == 1  
          initialdeltaE = deltaE;
          initialcumdE = cumdE;
          initialequicumdE = equicumdE;
          initialnewN = newN;
        end
      end

    % Apply smoothing of the path in CIELAB space if requested.  The aim is to
    % smooth out sharp lightness/colour changes that might induce the perception
    % of false features.  In doing this there will be some cost to the
    % perceptual contrast at these points.
    if sigma
      L = smooth(L, sigma, cyclic);
      a = smooth(a, sigma, cyclic);
      b = smooth(b, sigma, cyclic);
    end

    % Convert map back to RGB
    newlabmap = [L a b];
    newrgbmap = applycform(newlabmap, makecform('lab2srgb', 'AdaptedWhitePoint', wp));
    
    if diagnostics
        % Compute actual perceptual contrast values achieved
        if strcmpi(formula, 'CIE76')
          newdeltaE = cie76(L, a, b, W);
        elseif strcmpi(formula, 'CIEDE2000')
          newdeltaE = ciede2000(L, a, b, W);
        else
          error('Unknown colour difference formula')
        end

        show(sineramp,1,'Unequalised input colour map'), colormap(rgbmap);
        show(sineramp,2,'Equalised colour map'), colormap(newrgbmap);

        figure(3), clf
        subplot(2,1,1)
        plot(initialcumdE)
        title(sprintf('Cumulative change in %s of input colour map', formula));
        axis([1 N 0 1.05*max(cumdE(:))]);
        xlabel('Colour map index');

        subplot(2,1,2)
        plot(1:N, initialdeltaE, 1:N, newdeltaE)
        axis([1 N 0 1.05*max(initialdeltaE(:))]);
        legend('Original colour map', ...
         'Adjusted colour map', 'location', 'NorthWest');
        title(sprintf('Magnitude of %s differences along colour map', formula));
        xlabel('Colour map index');
        ylabel('dE')

        % Convert newmap back to Lab to check for gamut clipping
        labmap = applycform(newrgbmap,...
          makecform('srgb2lab', 'AdaptedWhitePoint', wp));
        figure(4), clf
        subplot(3,1,3)
        plot(1:N, L-labmap(:,1),...
         1:N, a-labmap(:,2),...
         1:N, b-labmap(:,3))
        legend('Lightness', 'a', 'b', 'Location', 'NorthWest')
        maxe = max(abs([L-labmap(:,1); a-labmap(:,2); b-labmap(:,3)]));
        axis([1 N -maxe maxe]);
        title('Difference between desired and achieved L a b values (gamut clipping)')
        
        % Plot RGB values
        subplot(3,1,1)
        plot(1:N, newrgbmap(:,1),'r-',...
         1:N, newrgbmap(:,2),'g-',...
         1:N, newrgbmap(:,3),'b-')

        legend('Red', 'Green', 'Blue', 'Location', 'NorthWest')
        axis([1 N 0 1.1]);        
        title('RGB values along colour map')
        
        % Plot Lab values
        subplot(3,1,2)
        plot(1:N, L, 1:N, a, 1:N, b)
        legend('Lightness', 'a', 'b', 'Location', 'NorthWest')
        axis([1 N -100 100]);                
        title('L a b values along colour map')

        
        %% Extra plots to generate figures.
        % The following plots were developed for illustrating the
        % equalization of MATLAB's hot(64) colour map.  Run using:
        % >> hoteq = equalisecolourmap('RGB', hot, 'CIE76', [1 0 0], 0, 0, 1);

        fs = 12; % fontsize
        % Colour difference values
        figure(10), clf
        plot(1:N, initialdeltaE, 'linewidth', 2);
        axis([1 N 0 1.05*max(initialdeltaE(:))]);
        title(sprintf('Magnitude of %s lightness differences along colour map', ...
          formula), 'fontweight', 'bold', 'fontsize', fs);
        
        xlabel('Colour map index', 'fontweight', 'bold', 'fontsize', fs);
        ylabel('dE', 'fontweight', 'bold', 'fontsize', fs);

        % Cumulative difference plot showing mapping of equicontrast colours
        figure(11), clf
        plot(initialcumdE, 'linewidth', 2), hold on
        s = 7;
        s = 10;
        for n = s:s:N
          line([0 initialnewN(n) initialnewN(n)],...
           [initialequicumdE(n) initialequicumdE(n) 0], ...
           'linestyle', '--', 'color', [0 0 0], 'linewidth', 1.5)

            ah = s/2; % Arrow height and width
            aw = ah/5;
            plot([initialnewN(n)-aw initialnewN(n) initialnewN(n)+aw],...
             [ah 0 ah],'k-', 'linewidth', 1.5)            
          end

          title(sprintf('Cumulative change in %s lightness of colour map', formula), ...
            'fontweight', 'bold', 'fontsize', fs);
          axis([1 N+6 0 1.05*max(cumdE(:))]);
          xlabel('Colour map index', 'fontweight', 'bold', 'fontsize', fs);
          ylabel('Equispaced cumulative contrast', 'fontweight', 'bold', 'fontsize', fs);
          hold off


        end

%----------------------------------------------------------------------------
%
% Function to smooth an array of values but also ensure end values are not
% altered or, if the map is cyclic, ensures smoothing is applied across the end
% points in a cyclic manner.  Assumed that the input data is a column vector

function Ls = smooth(L, sigma, cyclic)

  sze = ceil(6*sigma);  if ~mod(sze,2), sze = sze+1; end
  G = fspecial('gaussian', [sze 1], sigma);

  if cyclic 
        Le = [L(:); L(:); L(:)]; % Form a concatenation of 3 repetitions of the array. 
        Ls = filter2(G, Le);     % Apply smoothing filter
                                 % and then return the center section    
                                 Ls = Ls(length(L)+1: length(L)+length(L));

    else  % Non-cyclic colour map: Pad out input array L at both ends by 3*sigma
          % with additional values at the same slope.  The aim is to eliminate
          % edge effects in the filtering
          extension = (1:ceil(3*sigma))';

          dL1 = L(2)-L(1);
          dL2 = L(end)-L(end-1);
          Le = [-flipud(dL1*extension)+L(1); L;  dL2*extension+L(end)];

        Ls = filter2(G, Le);  % Apply smoothing filter
        
        % Trim off extensions
        Ls = Ls(length(extension)+1 : length(extension)+length(L));
      end

%----------------------------------------------------------------------------
% Delta E using the CIE76 formula + weighting

function deltaE = cie76(L, a, b, W)

  N = length(L);

    % Compute central differences 
    dL = zeros(size(L));
    da = zeros(size(a));
    db = zeros(size(b));

    dL(2:end-1) = (L(3:end) - L(1:end-2))/2;
    da(2:end-1) = (a(3:end) - a(1:end-2))/2;
    db(2:end-1) = (b(3:end) - b(1:end-2))/2;
    
    % Differences at end points
    dL(1) = L(2) - L(1);  dL(end) = L(end) - L(end-1);
    da(1) = a(2) - a(1);  da(end) = a(end) - a(end-1);
    db(1) = b(2) - b(1);  db(end) = b(end) - b(end-1);
    
    deltaE = sqrt(W(1)*dL.^2 + W(2)*da.^2 + W(3)*db.^2);
    
%----------------------------------------------------------------------------
% Delta E using the CIEDE2000 formula + weighting
%
% This function requires Gaurav Sharma's MATLAB implementation of the CIEDE2000
% color difference formula deltaE2000.m available at:
% http://www.ece.rochester.edu/~gsharma/ciede2000/

function deltaE = ciede2000(L, a, b, W)

  N = length(L);
  Lab = [L(:) a(:) b(:)];
  KLCH = 1./W;

    % Compute deltaE using central differences
    deltaE = zeros(N, 1);
    deltaE(2:end-1) = deltaE2000(Lab(1:end-2,:), Lab(3:end,:), KLCH)/2; 

    % Differences at end points    
    deltaE(1) = deltaE2000(Lab(1,:), Lab(2,:), KLCH);
    deltaE(end) = deltaE2000(Lab(end-1,:), Lab(end,:), KLCH);

% PBSPLINE - Basic Periodic B-spline
%
% Usage:  S = pbspline(P, k, N)
% 
% Arguments:   P - [dim x Npts] array of control points
%              k - order of spline (>= 2). 
%                  k = 2: Linear
%                  k = 3: Quadratic, etc
%              N - Optional number of points to evaluate along
%                  spline. Defaults to 100.
%
% Returns:     S - spline curve  [dim x N] spline points
%
% See also: BBSPLINE

% PK March 2014
%    Nov   2015  Made basis calculation slightly less wasteful

% Needs a bit of tidying up and checking on domain of curve.  Should be
% merged with BBSPLINE

function S = pbspline(P, k, N)

  if ~exist('N', 'var'), N = 100; end

    % For a closed spline check if 1st and last control points match.  If not
    % add another control point so that they do match
    if norm(P(:,1) - P(:,end)) > 0.01
      P = [P P(:,1)];
    end
    
    % Now add k - 1 control points that wrap over the first control points
    P = [P P(:,2:2+k-1)];
    
    [dim, np1] = size(P);
    n = np1-1;

    assert(k >= 2, 'Spline order must be 2 or greater');    
    assert(np1 >= k, 'No of control points must be >= k');
    assert(N >= 2, 'Spline must be evaluated at 2 or more points');
    
    % Form a uniform sequence. Number of knot points is m + 1 where m = n + k + 1
    ti = [0:(n+k+1)]/(n+k+1);
    nK = length(ti);

    % Domain of curve is [ti_k to ti_n] or [ti_(k+1) to ti_(n+1)] ???
    tstart = ti(k);
    tend = ti(n);

    dt = (tend-tstart)/(N-1);
    t = tstart:dt:tend;
    
    % Build complete array of basis functions.  We maintain two
    % arrays, one storing the basis functions at the current level of
    % recursion, and one storing the basis functions from the previous
    % level of recursion
    B = cell(1,nK-1);
    Blast = cell(1,nK-1);
    
    % 1st level of recursive construction
    for i = 1:nK-1
     Blast{i} = t >= ti(i) & t < ti(i+1)  & ti(i) < ti(i+1); 
   end

    % Subsequent levels of recursive basis construction.  Note the logic to
    % handle repeated knot values where ti(i) == ti(i+1)
    for ki = 2:k    
      for i = 1:nK-ki

        if (ti(i+ki-1) - ti(i)) < eps
          V1 = 0;
        else
          V1 = (t - ti(i))/(ti(i+ki-1) - ti(i)) .* Blast{i};
        end

        if (ti(i+ki) - ti(i+1)) < eps
          V2 = 0;
        else
          V2 = (ti(i+ki) - t)/(ti(i+ki) - ti(i+1)) .* Blast{i+1};
        end

        B{i} = V1 + V2;

%           This is the ideal equation that the code above implements  
%            N{i,ki} = (t - ti(i))/(ti(i+ki-1) - ti(i)) .* N{i,ki-1} + ...
%                      (ti(i+ki) - t)/(ti(i+ki) - ti(i+1)) .* N{i+1,ki-1};
end

        % Swap B and Blast, but only if this is not the last iteration
        if ki < k
          tmp = Blast;
          Blast = B;
          B = tmp;
        end
      end

    % Apply basis functions to the control points
    S = zeros(dim, length(t));
    for d = 1:dim
      for i = 1:np1
        S(d,:) = S(d,:) + P(d,i)*B{i};
      end
    end
    
   % Finally, because of the knot arrangements, the start of the spline may not
   % be close to the first control point if the spline order is 3 or greater.
   % Normally for a closed spline this is irrelevant.  However for our purpose
   % of using closed bplines to form paths in a colourspace this is important to
   % us.  The simple brute force solution used here is to search through the
   % spline points for the point that is closest to the 1st control point and
   % then rotate the spline points accordingly
   
   distsqrd = 0;
   for d = 1:dim
     distsqrd = distsqrd + (S(d,:) - P(d,1)).^2;
   end
   
   [~,ind] = min(distsqrd);
   
   S = circshift(S, [0, -ind+1]);

   function cmyk = rgb2cmyk(map)

    k = min(1-map,[],2);
    denom = 1 - k + eps;  % Add eps to avoid divide by 0
    c = (1-map(:,1) - k)./denom;
    m = (1-map(:,2) - k)./denom;
    y = (1-map(:,3) - k)./denom;
    
    cmyk = [c m y k];

    function Lab = rgb2lab(im, wp)

      if ~exist('wp', 'var'), wp = 'D65'; end

        cform = makecform('srgb2lab',...
          'adaptedwhitepoint', whitepoint(wp));    
        if strcmp(class(im),'uint8')
          im = double(im)/255;
        end
        Lab = applycform(im, cform);


        function nrgb = rgb2nrgb(im, offset)

          if ndims(im) ~= 3;
            error('Image must be a colour image');
          end

    % Convert to double if needed and define an offset = 1/255 max value to
    % be used in the normalization to avoid division by zero
    if ~strcmp(class(im), 'double')
      im = double(im);
      if ~exist('offset', 'var'), offset = 1; end
    else   % Assume we have doubles in range 0..1
      if ~exist('offset', 'var'), offset = 1/255; end
      end

      nrgb = zeros(size(im));
      gim = sum(im,3) + offset;

      nrgb(:,:,1) = im(:,:,1)./gim;
      nrgb(:,:,2) = im(:,:,2)./gim;
      nrgb(:,:,3) = im(:,:,3)./gim;            

      function n = normalise(im, reqmean, reqvar)

        if ~(nargin == 1 | nargin == 3)
         error('No of arguments must be 1 or 3');
       end

    if nargin == 1   % Normalise 0 - 1
  if ndims(im) == 3         % Assume colour image 
    hsv = rgb2hsv(im);
    v = hsv(:,:,3);
      v = v - min(v(:));    % Just normalise value component
      v = v/max(v(:));
      hsv(:,:,3) = v;
      n = hsv2rgb(hsv);
  else                      % Assume greyscale 
    if ~isa(im,'double'), im = double(im); end
      n = im - min(im(:));
      n = n/max(n(:));
    end

    else  % Normalise to desired mean and variance

  if ndims(im) == 3         % colour image?
    error('cannot normalise colour image to desired mean and variance');
  end

  if ~isa(im,'double'), im = double(im); end  
    im = im - mean(im(:));    
  im = im/std(im(:));      % Zero mean, unit std dev

  n = reqmean + im*sqrt(reqvar);
end

function [attributeStr,splineorder,colourspace,colpts] = map_sequence(rgb)
 % Function to map Andrew's RGB colorbrewer scheme
 attributeStr = 'linear';
 splineorder = 2;
 colourspace = 'RGB';
 c = rgb/255; 
 c = rgb2hsv(c);
 colpts = repmat(c,2,1);
 colpts(1,2:3) = [1 0.5];
 colpts(2,2:3) = [1,1];
 colpts = hsv2rgb(colpts);

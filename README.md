# Aquacrop_spatiallyexplicit
AquaCrop adaptations made for paper Van Loo &amp; Verstraeten 202X (submitted)

In this study , adaptions were made to the AquaCrop and run for the last 4000 years to 
simulate the impact of climate and land cover changes, as well as soil dynamics, 
on the productivity of winter wheat crops for a Mediterranean mountain environment in SW Turkey. 

# AquaCrop source code

The adaptations were based of AquaCropOS version 5.0a (Foster 2016).
Only files where adaptions were made are presented in this repository. For the complete AQuaCropOS code, please refer to their website.

# General code flow
```
├── START               
│
│
├── main_rout4os_WF.m              
│   ├── load in DTM    
│   │
|   ├── load in Slope    
│   │
|   ├── load in climate data    
│   │
|   ├── load in SWC regression    
│   │
|   ├── load in Lithology    
│   │
|   ├── load in Landuse    
│   │
│   └── load in Soil and Crop data structure
│
├── AQ_rout_mvl.m (or AQ_rout_mvl_spatialpar.m when using spatial parallelism)
│   │ 
│   ├── Initialize data matrices for output
│   │
|   ├──  Get runoff from neighbouring cells
│   │  
|   ├── Write input text files for AquaCrop run to disk  
│   │
│   └──Run AquaCropOS
│
└── AquaCropOS_RUN.m
    │ 
    ├── AquaCropOS calculations
    │
    ├──  adaptation of AOS_RainfallPartition.m
    │  
    └── AquaCropOS calculations (continued)
```
# Aquacrop adaptations
* (1) AquaCrop has been made spatially explicit allowing hydrological interactions between different
landscape positions. 
* (2) Runoff and re-infiltration processes were adapted to fit a
Mediterranean environment 
* (3) Computational time is kept limited by implementing parallelisation 
schemes on a supercomputer. 
  
 
## (1) AquaCrop spatially explicit
The main work is done inside `AQ_rout_mvl.m`. Here, cells in the catchment are ran in such order that whenever a given cell is calculated,
all upstream cells were already calculated. Hence, the potential runoff this given cell may receive is already calculated.
This information is stored in `extraPtemp`.

```Matlab
controlMrow=rowmatrix(rowY_select(j)-1:rowY_select(j)+1,colX_select(j)-1:colX_select(j)+1); % cel die je beschouwd, die zijn buren, naar waar geven ze hun afvoer? Y richting
controlMcol=colmatrix(rowY_select(j)-1:rowY_select(j)+1,colX_select(j)-1:colX_select(j)+1); % cel die je beschouwd, die zijn buren, naar waar geven ze hun afvoer? X richting

diffMrow=statusMrow-controlMrow;
diffMcol=statusMcol-controlMcol;

binaryM=zeros(size(statusMrow));
binaryM(diffMrow == 0 & diffMcol ==0)=1;

extraPM=extraPtemp(rowY_select(j)-1:rowY_select(j)+1,colX_select(j)-1:colX_select(j)+1,:); %what if this pixel is outside ctc? it shouldnt matter since its runoff will be zero and won't affect the sum of extraP (see RunTemp, outer pixels are empty)

binaryM=repmat(binaryM,[1,1,size(extraPM,3)]);
extraP4cell=binaryM.*extraPM;
extraP4cell=sum(sum(extraP4cell)); % 1e en 2e sum = x en y (2D), 3e dim blijft staan
extraP4cell=extraP4cell(:)'; %dagelijkse extra neerslag voor bepaald jaar, voor bepaalde cell (via runoff van omliggende pixels)
extraP4cellSum=extraP4cellSum+extraP4cell; %om op te tellen bij P voor de waterbalans => werkt dit in parfor?!
```
## (2) Runoff and re-infiltration

* `write_CN` calculates the CN based on the stoniness (relusting from erosion)
of the soil, Lithology (soil texture) and land cover
* `AOS_RainfallPartition.m` is adapted to take into account re-infiltration (see section 3.3.1 in paper)

## (3) Parallelisation schemes

### spatial-parallelism
When using spatial-parallelism, please refer to `AQ_rout_mvl_spatialpar.m` This replaces the original `AQ_rout_mvl.m` discussed untill now.

* Cells in the catchment are sorted ascending based on their upstream area.
* These cells are ran before all others with a larger upstream area 

```matlab
%% Preparation parfor loop
runoffcalcmasked=runoffcalc;
runoffcalcmasked([1 end],:)=0;
runoffcalcmasked(:,[1 end])=0;
[waarde,locatie]=sort(runoffcalcmasked(:),'ascend');
[rowY,colX]=ind2sub(size(runoffcalcmasked),locatie);
```
* a parforloop allows parralel computing inside MATLAB itself

`parfor j=1:it_num_amount
`
### time-parallelism

For this scheme, each year is ran in parallel. Equation 1 from the paper is used to 
estimate SWC at the start of each year. This is implemented in `main_rout4os.m`
 




# References 

1. Foster, T. (2016), AquaCropOS - Reference manual, Version 5.0a, Technical Report November, University of Manchester.
2. Schwanghart, W. and Kuhn, N. J. (2010), Environmental Modelling & Software TopoToolbox : A set of Matlab functions for topographic analysis, Environmental Modelling and Software 25(6), 770-781.
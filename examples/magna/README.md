# Overview

This is an example of associating data from the 2020 Magna, Utah 5.7 aftershock sequence.

## Step 1: Travel Time Tables

To begin we must generate travel time tables.  For simplicitly we will use a  1D model and the travel time calculator from [GrowClust](https://github.com/interseismic/PhaseLink/blob/master/raytracer.tar.gz).

The velocity model that we'll use in the Wasatch Front is given by depth (km), Vp (km/s), Vs (km/s) tuples

     0.0  3.40 1.95
    5.54  5.90 3.39
    21.1  6.40 3.68
    32.0  7.50 4.31
    46.0  7.90 4.54

Note, that this model has been modified so that 0 is actually +2km above sea-level.  Call this file wasatch.txt.

Next, we run GrowClust's travel time calculator with the following script

    #!/bin/bash
    ttimes_growclust=make_TTable
    velocity_model="wasatch.txt"
    p_ttable="TT.wasatch.pg"
    s_ttable="TT.wasatch.sg"
    pphase=1
    sphase=2
    kmsec=1
    # Table depth spacing: Min/Max/Increment
    depth_spacing="0.0 40.0 0.5"
    # Table distance spacing: Min/Max/Increment
    distance_spacing="0 200 1"
    # Background vp/vs
    vpvs=1.74
    # Ray parameter bound - we use 32 km as the mangle which corresponds
    # to 1/7.5 and 1/4.31 as seen in the following 2 lines 
    noPn=0.133
    noSn=0.238
    # First column is depth
    vmodel_first_column_format=1
    # Read the velocity model into a
    n_lines_in_model=5
    # Read the file
    vmodel=`cat ${velocity_model}`
    n_lines_in_model=`wc -l ${velocity_model}`
    echo "$vmodel"

    ${ttimes_growclust} << EOF
    ${p_ttable}
    ${pphase}
    ${kmsec}
    ${depth_spacing}
    ${distance_spacing}
    ${vpvs}
    ${noPn}
    ${velocity_model}
    ${vmodel_first_column_format}
    ${n_lines_in_model}
    ${vmodel}
    EOF

    ${ttimes_growclust} << EOF
    ${s_ttable}
    ${sphase}
    ${kmsec}
    ${depth_spacing}
    ${distance_spacing}
    ${vpvs}
    ${noSn}
    ${velocity_model}
    ${vmodel_first_column_format}
    ${n_lines_in_model}
    ${vmodel}
    EOF

After creating the travel time tables we will repackage them into a 3D point cloud for use by the associator.  To do this run

    ttables


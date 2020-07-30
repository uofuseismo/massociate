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

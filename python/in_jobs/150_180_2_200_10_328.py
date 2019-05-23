#!/usr/bin/env python3
#
# Incoming propagation job: v1.0
# Time of writing: 2019-05-22 17:30:49.645972
#
# Platform
# Node=Hyper-Abacus
# Machine=AMD64
# System=Windows
# Version=10.0.18362
# Release=10
# Processor=Intel64 Family 6 Model 69 Stepping 1, GenuineIntel
#
# Setup
# Runs=2
# Seed=None
# Outgoing_File=./out_jobs/telemetry/150_180_2_200_10_328.outgoing
#
# Parameters
# Z=2 [proton number]
# A=None [atomic mass units]
# E=2e+20 [electron volts]
# Algorithm=dop853
# Max_Step=0.01 [AU]
# R_Limit=6.0 [AU]
# B_Override=None [T]
# Step_Override=None [AU]
# Origin=[ 9.99978683e-01  2.61063868e-21 -3.69229629e-05] [AU]
# Position=[-1.21932665  0.61683717 -5.84234775] [AU]
# Beta=[ 0.35339257 -0.09822249  0.93030427]
#
# Script

import gz

incoming = gz.path.Incoming([0.9999786825174377, 2.6106386785763696e-21, -3.692296288731554e-05], [-1.2193266496489956, 0.6168371718512692, -5.842347745406368], [0.35339257073154473, -0.09822248534202761, 0.9303042697553233], 2, None, 2e+20, 6.174316738356785, max_step=0.01, R_limit=6.0, save_path=None, filename='150_180_2_200_10_328_0000')
incoming.propagate(algorithm='dop853')

incoming = gz.path.Incoming([0.9999786825174377, 2.6106386785763696e-21, -3.692296288731554e-05], [-1.2193266496489956, 0.6168371718512692, -5.842347745406368], [0.35339257073154473, -0.09822248534202761, 0.9303042697553233], 2, None, 2e+20, 5.154166302928779, max_step=0.01, R_limit=6.0, save_path=None, filename='150_180_2_200_10_328_0001')
incoming.propagate(algorithm='dop853')


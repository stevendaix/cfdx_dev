#!/usr/bin/env python3
# Code_Saturne Python case definition
case.set('physics_model', 'navier_stokes')
case.set('turbulence_model', 'KEPSILON')
case.set('energy_model', True)
case.set('transient', False)
case.set('mach_number', 0.8)
case.set('alpha', 1.25)
case.set('freestream_pressure', 101325.0)
case.set('freestream_temperature', 288.15)
case.set('material', 'air')
case.set('boundary_inlet_type', 'inlet')
case.set('boundary_outlet_type', 'outlet')
case.set('boundary_wall_type', 'wall')

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 09:50:37 2023

@author: kutsj
"""
import numpy as np
import pyvtk

import SpatialFiber as sf

def  fibrous_system(properties, outputf_file):
    pass
    



def main():
    
    rve_prop = {
        "num_fibers": 100,
        "radius": 1,
        "height": 50,
        "width": 50,
        "depth": 10,
        "center": (0,0,0),
        "tolerance": 0
    }
    
    fibrous_system(rve_prop, "rve1.vtk")

    
if __name__=="__main__":
    main()


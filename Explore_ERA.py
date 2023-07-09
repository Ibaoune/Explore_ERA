#!/usr/bin/env python3
# -*- coding: utf-8 -*-



__doc__ = """
Created on Wed Jun 21 14:57:44 2023
"""
__author__ = "M. El Aabbaribaoune (@_um6p)"

import os
import sys
import json
import core

def main(pltcfgfile):
    
    # Load the configuration
    with open(pltcfgfile) as file:
        plotcfg = json.load(file)
        
    inputdir = plotcfg['inputdir']
    
    if not os.path.exists(inputdir):
        print(inputdir," does not exists. Exit.")
        sys.exit()
        
    os.chdir(inputdir)
    print("Working inside ",inputdir)
    
            
    if "plotWhenItRains" in plotcfg:
            
        core.plot_flds_when_it_rains(plotcfg, inputdir)
           
            
if __name__ == "__main__":
    
    pltcfgfile = '/home/mohammad/Desktop/Mohammad/Pro/PostPhD/UM6P/Y_2023/others/for_Mme_Ouaraini/08-07-2023/explore_era.json'
    # Read cfg file from user input if given. It overrides default
    #if len(sys.argv) > 1: pltcfgfile = sys.argv[1]
    
    main(pltcfgfile)  
    
    

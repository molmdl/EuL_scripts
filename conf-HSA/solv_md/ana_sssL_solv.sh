#!/bin/bash

python ../ana_code/metal_geo_analysis.py --system phe_sssL_sap --tpr phe_sssL_sap/prod_0.tpr --xtc phe_sssL_sap/prod_4x100ns.xtc --outdir merged 
python ../ana_code/metal_geo_analysis.py --system phe_sssL_sap --tpr phe_sssL_sap/prod_1.tpr --xtc phe_sssL_sap/prod_1.xtc --outdir t0 
python ../ana_code/metal_geo_analysis.py --system phe_sssL_sap --tpr phe_sssL_sap/prod_1.tpr --xtc phe_sssL_sap/prod_1.xtc --outdir t1 
python ../ana_code/metal_geo_analysis.py --system phe_sssL_sap --tpr phe_sssL_sap/prod_2.tpr --xtc phe_sssL_sap/prod_2.xtc --outdir t2 
python ../ana_code/metal_geo_analysis.py --system phe_sssL_sap --tpr phe_sssL_sap/prod_3.tpr --xtc phe_sssL_sap/prod_3.xtc --outdir t3 
python ../ana_code/metal_geo_analysis.py --system me_sssL_sap --tpr me_sssL_sap/prod_0.tpr --xtc me_sssL_sap/prod_4x100ns.xtc --outdir merged 
python ../ana_code/metal_geo_analysis.py --system me_sssL_sap --tpr me_sssL_sap/prod_0.tpr --xtc me_sssL_sap/prod_0.xtc --outdir t0 
python ../ana_code/metal_geo_analysis.py --system me_sssL_sap --tpr me_sssL_sap/prod_1.tpr --xtc me_sssL_sap/prod_1.xtc --outdir t1 
python ../ana_code/metal_geo_analysis.py --system me_sssL_sap --tpr me_sssL_sap/prod_2.tpr --xtc me_sssL_sap/prod_2.xtc --outdir t2 
python ../ana_code/metal_geo_analysis.py --system me_sssL_sap --tpr me_sssL_sap/prod_3.tpr --xtc me_sssL_sap/prod_3.xtc --outdir t3

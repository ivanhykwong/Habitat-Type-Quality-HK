import arcpy
from arcpy import env
from arcpy.sa import *
import os
import glob

# copy raster

os.chdir("D:/landsat/raw") 
arcpy.env.overwriteOutput = True
env.workspace = "D:/landsat/raw"
tif_list = glob.glob("*.tif")

for t in tif_list:
    arcpy.management.CopyRaster(t, "D:/landsat/preprocess/" + t, '', None, "-3.402823e+38", "NONE", "NONE", '32_BIT_FLOAT', "NONE", "NONE", "TIFF", "NONE", "CURRENT_SLICE", "NO_TRANSPOSE")


# Move problematic scenes to another folder

import shutil
os.chdir("D:/landsat/preprocess/c")
with open("D:/landsat/scenewithproblem.txt") as file: 
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]

lines = [l[:-4] for l in lines]
tif_list = glob.glob("*")

for t in lines:
    t_list = [s for s in tif_list if t in s] # filter if file name contain target
    for f in t_list:
        shutil.move(f, "problem/")


# Find unique date & merge

os.chdir("D:/landsat/preprocess/")
arcpy.env.overwriteOutput = True
env.workspace = "D:/landsat/preprocess/"
tif_list = glob.glob("L*.tif")

from collections import OrderedDict 
tif_date = [t[17:25] for t in tif_list]
col_date = list(OrderedDict.fromkeys(tif_date)) # find unique
for d in col_date:
    d_list = [s for s in tif_list if d in s] # filter date if contain
    if len(d_list) == 1:
        arcpy.management.CopyRaster(d_list[0], "D:/landsat/preprocess/m/" + d_list[0], '', None, "-3.402823e+38", "NONE", "NONE", '', "NONE", "NONE", "TIFF", "NONE", "CURRENT_SLICE", "NO_TRANSPOSE")
    else:
        # change no of bands MSS 4  TM 6  ETM 6  OLI 7
        arcpy.management.MosaicToNewRaster(d_list, "m", d_list[0], None, "32_BIT_FLOAT", 30, 6, "MEAN", "FIRST") 


# Reproject to HK1980

os.chdir("D:/landsat/preprocess/m")
arcpy.env.overwriteOutput = True
env.workspace = "D:/landsat/preprocess/m"
tif_list = glob.glob("*.tif")
for t in tif_list:
    arcpy.management.ProjectRaster(t, t[:-4]+"_hk1980.tif", 'PROJCS["Hong_Kong_1980_Grid",GEOGCS["GCS_Hong_Kong_1980",DATUM["D_Hong_Kong_1980",SPHEROID["International_1924",6378388.0,297.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",836694.05],PARAMETER["False_Northing",819069.8],PARAMETER["Central_Meridian",114.1785555555556],PARAMETER["Scale_Factor",1.0],PARAMETER["Latitude_Of_Origin",22.31213333333334],UNIT["Meter",1.0]]', "BILINEAR", "30 30", "Hong_Kong_1980_To_WGS_1984_1", "798390 799500")

# Clip rectangle AOI

os.chdir("D:/landsat/preprocess/hk1980")
arcpy.env.overwriteOutput = True
env.workspace = "D:/landsat/preprocess/hk1980"
tif_list = glob.glob("*.tif")

for t in tif_list:
    arcpy.management.Clip(t, "798390 799500 865740 848940", "D:/landsat/preprocess/rectangle/" + t, None, "-3.402823e+38", "NONE", "MAINTAIN_EXTENT")


# mask cloud cover for MSS

os.chdir("D:/landsat/preprocess/MSS") 
arcpy.env.overwriteOutput = True
env.workspace = "D:/landsat/preprocess/MSS"
tif_list = glob.glob("*.tif")
for t in tif_list:
    shp_name = "D:/landsat/preprocess/mss_cloud/"+t[:-4]+".shp"
    shp_erase = "D:/landsat/preprocess/MSS_aoi.shp"
    if os.path.isfile(shp_name) == True:
        arcpy.analysis.Erase("D:/landsat/aoi/rectangle.shp", shp_name, shp_erase, None)
        arcpy.management.Clip(t, "798390 799500 865740 848940", "D:/landsat/preprocess/MSS/cloudmask/" + t, shp_erase, "-3.402823e+38", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    else:
        arcpy.management.CopyRaster(t, "D:/landsat/preprocess/MSS/cloudmask/" + t, '', None, "-3.402823e+38", "NONE", "NONE", '', "NONE", "NONE", "TIFF", "NONE", "CURRENT_SLICE", "NO_TRANSPOSE")



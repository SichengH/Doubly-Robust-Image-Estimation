#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 11:48:50 2019

@author: sichenghao
"""
import os
import pylidc as pl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manim
from skimage.measure import find_contours
import pylidc as pl
import pdb
from pylidc.utils import consensus
from PIL import Image
os.getcwd()
#%%
# Query for all CT scans with desired traits.
#scans = pl.query(pl.Scan).filter(pl.Scan.slice_thickness <= 1,pl.Scan.pixel_spacing <= 0.6)
scans = pl.query(pl.Scan)
print(scans.count())
#%%
def make_image(data, outputname, size=(10, 10), dpi=80):
    fig = plt.figure()
    fig.set_size_inches(size)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.set_cmap('gray')
    ax.imshow(data, aspect='equal')
    plt.savefig(outputname, dpi=dpi)

def crop_center(img,cropx,cropy):
    y,x = img.shape 
    startx = x//2-(cropx//2) 
    starty = y//2-(cropy//2)
    return img[starty:starty+cropy,startx:startx+cropx]
#%%
for i in range(2):
    scan = scans[i]
    nods = scan.cluster_annotations()
    #print(len(nods))
    print(i)
    count = 0
    for j in range(len(nods)):
        anno = nods[j]
        for k in range(len(anno)):
            ann = scan.annotations[count]
            count+=1
            #print(ann.Malignancy)
            #save mask only 
            padding = [(30,30), (30,30), (0,0)]
            mask = ann.boolean_mask(pad=padding)
            bbox = ann.bbox(pad=padding)
            img = mask[:,:,0]
            crop_size = min(img.shape[:2])
            img_crop = crop_center(img,crop_size,crop_size).astype(np.float32)
            gray_img = (((img_crop - img_crop.min()) / (img_crop.max() - img_crop.min())) * 255.9).astype(np.uint8) 
            gray_img = Image.fromarray(gray_img) 
            
            
#            fig,ax = plt.subplots(1,figsiz1e=(5,3))
#            plt.imshow(img_crop, cmap=plt.cm.gray)
#            plt.tight_layout()
#            plt.axis('off')
            IMG_name = str(ann.scan)
            pid = IMG_name.split(",")[1]
            pid = pid.split("=")[1]
            pid = pid.split(")")[0]
            pid = pid.split("-")[2]
            iid = pid+"_"+str(j)+"_"+str(ann.id)
            IMG_PATH = os.path.join("/Users/sichenghao/Desktop/DCM-Processing/images",iid)
            gray_img.save(IMG_PATH+'.png')
            #print(IMG_PATH)
 #           plt.tight_layout()
#            plt.savefig(IMG_PATH, dpi=80)
            
#            pdb.set_trace()
#            size = img_crop.shape[:2]
#            make_image(img_crop,IMG_PATH,size)
#            print(size)
#            plt.show()
#%%





pid = 'LIDC-IDRI-0001'
scan = pl.query(pl.Scan).filter(pl.Scan.patient_id == pid).first()
nods = scan.cluster_annotations()
print("%s has %d nodules." % (scan, len(nods)))
# => Scan(id=1,patient_id=LIDC-IDRI-0078) has 4 nodules.
for i,nod in enumerate(nods):
    print("Nodule %d has %d annotations." % (i+1, len(nods[i])))
# => Nodule 1 has 4 annotations.
# => Nodule 2 has 4 annotations.
# => Nodule 3 has 1 annotations.
# => Nodule 4 has 4 annotations.
    
vol = scan.to_volume()
print(vol.shape)
# => (512, 512, 87)

print("%.2f, %.2f" % (vol.mean(), vol.std()))
# => -702.15, 812.52

scan.visualize(annotation_groups=nods)
#

ann = pl.query(pl.Annotation).first()
ann = scan.annotations[0]
                                   
print(ann.scan.patient_id)

anns = pl.query(pl.Annotation).filter(pl.Annotation.spiculation == 5,
                                      pl.Annotation.malignancy == 5)
print(anns.count())

print(type(pl.Annotation.lobulation))
# => <class 'sqlalchemy.orm.attributes.InstrumentedAttribute'>
# => ^^^ queryable because it's an sqlalchemy attribute.

# Whereas ...
print(type(pl.Annotation.Lobulation))
# => <type 'property'>
# ^^^ not queryable because it's a computed property

ann = pl.query(pl.Annotation)\
        .filter(pl.Annotation.malignancy == 5).first()

print(ann.malignancy, ann.Malignancy)
# => 5, 'Highly Suspicious'

print(ann.margin, ann.Margin)
# => 2, 'Near Poorly Defined'
#('subtlety',
# 'internalStructure',
#'calcification',
# 'sphericity',
# 'margin',
# 'lobulation',
# 'spiculation',
# 'texture',
# 'malignancy')

ann.print_formatted_feature_table()

 #Query
svals = pl.query(pl.Annotation.spiculation)\
          .filter(pl.Annotation.spiculation > 3)

print(svals[0])
# => (4,)

print(all([s[0] > 3 for s in svals]))


ann = pl.query(pl.Annotation).first()
contours = ann.contours

print(contours[0])
# => Contour(id=21,annotation_id=1)

print("%.2f mm, %.2f mm^2, %.2f mm^3" % (ann.diameter,
                                         ann.surface_area,
                                         ann.volume))
# => 20.84 mm, 1242.74 mm^2, 2439.30 mm^3


mask = ann.boolean_mask()
print(mask.shape, mask.dtype)
# => (34, 27, 6), dtype('bool')

bbox = ann.bbox()
print(bbox)
# => (slice(151, 185, None), slice(349, 376, None), slice(44, 50, None))

vol = ann.scan.to_volume()
print(vol[bbox].shape)
# => (34, 27, 6)

print(ann.bbox_dims())
# => [21.45, 16.90, 15.0]




#vis



ann = pl.query(pl.Annotation).first()
vol = ann.scan.to_volume()

padding = [(30,30), (30,30), (0,0)]

mask = ann.boolean_mask(pad=padding)
bbox = ann.bbox(pad=padding)

fig,ax = plt.subplots(1,2,figsize=(5,3))

ax[0].imshow(vol[bbox][:,:,4], cmap=plt.cm.gray)
ax[0].axis('on')

ax[1].imshow(mask[:,:,2], cmap=plt.cm.gray)
ax[1].axis('on')

plt.tight_layout()
plt.savefig("./images/mask_bbox.png", bbox_inches="tight")
plt.show()



#save mask only 
padding = [(30,30), (30,30), (0,0)]
mask = ann.boolean_mask(pad=padding)
bbox = ann.bbox(pad=padding)
fig,ax = plt.subplots(1,figsize=(5,3))
ax.imshow(mask[:,:,2], cmap=plt.cm.gray)
ax.axis('off')

IMG_name = str(ann.scan)
sid = IMG_name.split(",")[0]
pid = IMG_name.split(",")[1]
sid = sid.split("=")[1]
pid = pid.split("=")[1]
pid = pid.split(")")[0]
pid = pid.split("-")[2]
iid = sid + "_"+pid
IMG_PATH = os.path.join("./images",iid)
plt.tight_layout()
plt.savefig(IMG_PATH, bbox_inches="tight")
plt.show()

ann = pl.query(pl.Annotation)\
        .filter(pl.Annotation.lobulation == 5).first()
ann.visualize_in_3d()

ann = pl.query(pl.Annotation).first()
ann = scan.annotations[0]
ann.visualize_in_scan()





# Query for a scan, and convert it to an array volume.
scan = pl.query(pl.Scan).filter(pl.Scan.patient_id == 'LIDC-IDRI-0078').first()
vol = scan.to_volume()

# Cluster the annotations for the scan, and grab one.
nods = scan.cluster_annotations()
anns = nods[0]

# Perform a consensus consolidation and 50% agreement level.
# We pad the slices to add context for viewing.
cmask,cbbox,masks = consensus(anns, clevel=0.5,
                              pad=[(20,20), (20,20), (0,0)])

# Get the central slice of the computed bounding box.
k = int(0.5*(cbbox[2].stop - cbbox[2].start))

# Set up the plot.
fig,ax = plt.subplots(1,1,figsize=(5,5))
ax.imshow(vol[cbbox][:,:,k], cmap=plt.cm.gray, alpha=0.5)

# Plot the annotation contours for the kth slice.
colors = ['r', 'g', 'b', 'y']
for j in range(len(masks)):
    for c in find_contours(masks[j][:,:,k].astype(float), 0.5):
        label = "Annotation %d" % (j+1)
        plt.plot(c[:,1], c[:,0], colors[j], label=label)

# Plot the 50% consensus contour for the kth slice.
for c in find_contours(cmask[:,:,k].astype(float), 0.5):
    plt.plot(c[:,1], c[:,0], '--k', label='50% Consensus')

ax.axis('off')
ax.legend()
plt.tight_layout()
#plt.savefig("../images/consensus.png", bbox_inches="tight")
plt.show()






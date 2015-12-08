#take numpy array; save an image with the indicated color table.
#based on arraytoimage singleColorImage, but adjusted to avoid throwing an error when a colortable is used

# Creates a false-color image of the array data.
# The default color table is gray.

from PIL import Image

import numpy

def saveimg(data, filename, colortable = None, autoscale = True, xsize = None, ysize = None):

    # default grayscale color table
    r = numpy.arange(0, 256, dtype = numpy.uint8)
    g = numpy.arange(0, 256, dtype = numpy.uint8)
    b = numpy.arange(0, 256, dtype = numpy.uint8)
    if colortable.any():#This line is changed from Hazen's code. otherwise it is the same.
        r = colortable[0].astype(numpy.uint8)
        g = colortable[1].astype(numpy.uint8)
        b = colortable[2].astype(numpy.uint8)

    # create the image
    temp = data.copy().astype(numpy.float)
    if autoscale:
        temp = temp - numpy.min(temp)
        if(numpy.max(temp)>0):
            temp = 255.0 * temp/numpy.max(temp)
    temp = numpy.round(temp)
    img = numpy.zeros((data.shape[0], data.shape[1], 4), numpy.uint8)
    img[:,:,0] = r[temp.astype(numpy.uint8)]
    img[:,:,1] = g[temp.astype(numpy.uint8)]
    img[:,:,2] = b[temp.astype(numpy.uint8)]

    # this sets the image transparency to zero
    img[:,:,3] = 255

    pilimage = Image.fromarray(img, "RGBA")

    if xsize and ysize:
        pilimage = pilimage.resize([xsize, ysize])
    pilimage.save(filename + ".png")
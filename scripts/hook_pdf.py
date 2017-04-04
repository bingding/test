#!/usr/bin/python
from __future__ import division

import numpy as np
import popplerqt4

def convertQImageToMat(incomingImage):
    '''  Converts a QImage into an opencv MAT format  '''

    incomingImage = incomingImage.convertToFormat(4)

    width = incomingImage.width()
    height = incomingImage.height()

    ptr = incomingImage.bits()
    ptr.setsize(incomingImage.byteCount())
    arr = np.array(ptr).reshape(height, width, 4)  #  Copies the data
    return arr


def load_pdf_file(pdf_file, res=144):
    """Arrays of pixels for each page with a given resolution res"""
    d = popplerqt4.Poppler.Document.load(pdf_file)
    pages = []
    for page in range(d.numPages()):
        qimage = d.page(page).renderToImage(res, res)
        pages.append(convertQImageToMat(qimage))
        
    return pages

def compare_pdfs(file1, file2):
    doc1 = load_pdf_file(file1)
    doc2 = load_pdf_file(file2)
    # Break if not same number of pages
    if len(doc1) != len(doc2):
        return False
    diff = [np.abs(np.sum(p1 - p2)) / p1.size for p1, p2 in zip(doc1, doc2)]
    if np.sum(diff) / len(p1) <= 0.:
        return True
    
    print "PDF mismatch: ",
    for kk, dd in enumerate(diff):
        print "page {}: {:.3}% ".format(kk, dd*100),
    print 
    return False

if __name__ == "__main__":
    pages = load_pdf_file("A4180-DV16-7-P19-BC89-DS_S1_L1234_160912_130536_cn.pdf")
    print pages[0].dtype
    

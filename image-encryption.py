#### IMPORTS ####

from PIL import Image
import os
import tkinter as tk
from tkinter import filedialog
import hashlib
import binascii
import textwrap
import cv2
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from importlib import reload
from bisect import bisect_left as bsearch
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ThreadPoolExecutor, as_completed


dna={}
dna["00"]="A"
dna["01"]="T"
dna["10"]="G"
dna["11"]="C"
dna["A"]=[0,0]
dna["T"]=[0,1]
dna["G"]=[1,0]
dna["C"]=[1,1]
#DNA xor
dna["AA"]=dna["TT"]=dna["GG"]=dna["CC"]="A"
dna["AG"]=dna["GA"]=dna["TC"]=dna["CT"]="G"
dna["AC"]=dna["CA"]=dna["GT"]=dna["TG"]="C"
dna["AT"]=dna["TA"]=dna["CG"]=dna["GC"]="T"
# Maximum time point and total number of time points
tmax, N = 100, 10000


# Function to select an image from the local file system
def image_selector():
    Tk().withdraw()  # Prevents the root window from appearing
    image_path = askopenfilename(title="Select an image", filetypes=[("Image files", "*.jpeg;*.jpg;*.png")])
    
    if image_path:
        image_type = imghdr.what(image_path)  # Check the image type
        if image_type in ['jpeg', 'jpg', 'png']:  # Check if the image is JPEG/JPG or PNG
            print("Image loaded:", image_path)
            return image_path
        else:
            print("Error: Unsupported image format! Please upload a JPEG/JPG or PNG image.")
            return None
    else:
        print("Error: No image selected!")
        return None

# Function to split an image into RGB channels
def split_into_rgb_channels(image):
    red = image[:, :, 2]
    green = image[:, :, 1]
    blue = image[:, :, 0]
    return red, green, blue

def securekey(iname):
    # Check if the file exists before attempting to open it
    if not os.path.isfile(iname):
        return None, None, None  # Return None values if the file doesn't exist

    img = Image.open(iname)
    m, n = img.size
    print("pixels: {0}  width: {2} height: {1} ".format(m * n, m, n))
    pix = img.load()
    plainimage = bytearray()

    for y in range(n):
        for x in range(m):
            pixel = pix[x, y]
            if isinstance(pixel, int):  # Handle grayscale or non-RGB images
                plainimage.extend([pixel, pixel, pixel])
            else:
                for value in pixel:
                    plainimage.append(value)

    key = hashlib.sha256()
    key.update(plainimage)
    return key.hexdigest(), m, n

x = 0.1
y = 0.1
z = 0.1

def baker_map(x, y, z):
  a = 1000
  b = 1000
  x_next = (2 * x) % a
  y_next = (y + x_next) % b
  z_next = z
  return x_next, y_next, z_next


def generate_bakers_map(key, width, height):
    x = 0.1
    y = 0.1
    z = 0.1
    # Generate sequences for 'x', 'y', and 'z'
    x_sequence = []
    y_sequence = []
    z_sequence = []

    for i in range(width * height):
        x, y, z = baker_map(x, y, z)
        x_sequence.append(x)
        y_sequence.append(y)
        z_sequence.append(z)

    return x_sequence, y_sequence, z_sequence

# Assuming securekey and selected_image_path functions are defined elsewhere
key, width, height = securekey(selected_image_path)
x_sequence, y_sequence, z_sequence = generate_bakers_map(key,width, height)

def update_baker_2d(key, x, y, z):
    x = 0.1
    y = 0.1
    z = 0.1
    key_bin = bin(int(key, 16))[2:].zfill(256)
    key_32_parts = [key_bin[i:i + 8] for i in range(0, len(key_bin), 8)]
    k = np.array([int(part, 2) for part in key_32_parts])

    t1 = k[:11].sum() ^ k[11:22].sum() ^ k[22:32].sum()

    x_next = 2 * x if x < 0.5 else 2 * x - 1
    y_next = 0.5 * y if x < 0.5 else 0.5 * y + 0.5
    z_next = z + t1 / 256

    return x_next, y_next, z_next


def decompose_matrix(iname):
    image = cv2.imread(iname)
    blue,green,red = split_into_rgb_channels(image)
    for values, channel in zip((red, green, blue), (2,1,0)):
        img = np.zeros((values.shape[0], values.shape[1]), dtype = np.uint8)
        img[:,:] = (values)
        if channel == 0:
            B = np.asmatrix(img)
        elif channel == 1:
            G = np.asmatrix(img)
        else:
            R = np.asmatrix(img)
    return B,G,R

def dna_encode(b, g, r):
    b = np.unpackbits(b, axis=1)
    g = np.unpackbits(g, axis=1)
    r = np.unpackbits(r, axis=1)

    m, n = b.shape
    enc_shape = (m, n // 2)

    def encode_channel(channel):
        encoded = np.empty(enc_shape, dtype='U1')
        for j in range(m):
            encoded[j] = np.array([dna[f"{channel[j, i]}{channel[j, i+1]}"] for i in range(0, n, 2)])
        return encoded

    b_enc = encode_channel(b)
    g_enc = encode_channel(g)
    r_enc = encode_channel(r)

    return b_enc, g_enc, r_enc

def key_matrix_encode(key,b):
    #encoded key matrix
    b = np.unpackbits(b,axis=1)
    m,n = b.shape
    key_bin = bin(int(key, 16))[2:].zfill(256)
    Mk = np.zeros((m,n),dtype=np.uint8)
    x=0
    for j in range(0,m):
            for i in range(0,n):
                Mk[j,i]=key_bin[x%256]
                x+=1

    Mk_enc=np.chararray((m,int(n/2)))
    idx=0
    for j in range(0,m):
        for i in range(0,n,2):
            if idx==(n/2):
                idx=0
            Mk_enc[j,idx]=dna["{0}{1}".format(Mk[j,i],Mk[j,i+1])]
            idx+=1
    Mk_enc=Mk_enc.astype(str)
    return Mk_enc

def xor_operation(b, g, r, mk):
    m, n = b.shape

    b_mk = np.array([dna[f"{b[i, j]}{mk[i, j]}"] for i in range(m) for j in range(n)])
    g_mk = np.array([dna[f"{g[i, j]}{mk[i, j]}"] for i in range(m) for j in range(n)])
    r_mk = np.array([dna[f"{r[i, j]}{mk[i, j]}"] for i in range(m) for j in range(n)])

    bx = np.chararray((m, n))
    gx = np.chararray((m, n))
    rx = np.chararray((m, n))

    bx[:, :] = b_mk.reshape((m, n))
    gx[:, :] = g_mk.reshape((m, n))
    rx[:, :] = r_mk.reshape((m, n))

    return bx, gx, rx


def gen_chaos_seq(m, n):
    global x0, y0, z0
    N = m * n * 4

    key_bin = bin(int(key, 16))[2:].zfill(256)
    key_32_parts = np.array([key_bin[i:i+8] for i in range(0, len(key_bin), 8)], dtype=np.uint8)
    k = np.packbits(key_32_parts, axis=0)

    k_sum = np.zeros(32, dtype=np.uint8)
    for i in range(32):
        k_sum[i] = np.packbits(k[i::32]).sum()

    k_sum_all = np.tile(k_sum, N // 32 + 1)[:N]

    x, y, z = np.zeros(N), np.zeros(N), np.zeros(N)

    for i in range(N):
        x0 = x0 + k_sum_all[i] / 256
        y0 = y0 + k_sum_all[i] / 256
        z0 = z0 + k_sum_all[i] / 256
        x[i] = x0
        y[i] = y0
        z[i] = z0
    return x, y, z

def parallel_sequence_indexing(x, y, z):
    n = len(x)

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(np.argsort, arr) for arr in [x, y, z]]

    k2_x = np.empty(n, dtype=np.uint32)
    k2_y = np.empty(n, dtype=np.uint32)
    k2_z = np.empty(n, dtype=np.uint32)

    results = list(as_completed(futures))

    k2_x[results[0].result()] = np.arange(n, dtype=np.uint32)
    k2_y[results[1].result()] = np.arange(n, dtype=np.uint32)
    k2_z[results[2].result()] = np.arange(n, dtype=np.uint32)

    return k2_x, k2_y, k2_z


def scramble(fx, fy, fz, b, g, r):
    p, q = b.shape
    size = p * q

    bx = b.reshape(size).astype(str)
    gx = g.reshape(size).astype(str)
    rx = r.reshape(size).astype(str)

    bx_s = np.empty(size, dtype='U1')
    gx_s = np.empty(size, dtype='U1')
    rx_s = np.empty(size, dtype='U1')

    for i in range(size):
        idx = fz[i]
        bx_s[i] = bx[idx]
    for i in range(size):
        idx = fy[i]
        gx_s[i] = gx[idx]
    for i in range(size):
        idx = fx[i]
        rx_s[i] = rx[idx]

    b_s = bx_s.reshape(p, q)
    g_s = gx_s.reshape(p, q)
    r_s = rx_s.reshape(p, q)
    return b_s, g_s, r_s


def decryption_scrambling(fx, fy, fz, b, g, r):
    p, q = b.shape
    size = p * q

    bx = b.flatten()
    gx = g.flatten()
    rx = r.flatten()

    bx_s = np.empty(size, dtype='U1')
    gx_s = np.empty(size, dtype='U1')
    rx_s = np.empty(size, dtype='U1')

    bx_s[fz] = bx
    gx_s[fy] = gx
    rx_s[fx] = rx

    b_s = bx_s.reshape((p, q))
    g_s = gx_s.reshape((p, q))
    r_s = rx_s.reshape((p, q))

    return b_s, g_s, r_s

def dna_decode(b,g,r):
    m,n = b.shape
    r_dec= np.ndarray((m,int(n*2)),dtype=np.uint8)
    g_dec= np.ndarray((m,int(n*2)),dtype=np.uint8)
    b_dec= np.ndarray((m,int(n*2)),dtype=np.uint8)
    for color,dec in zip((b,g,r),(b_dec,g_dec,r_dec)):
        for j in range(0,m):
            for i in range(0,n):
                dec[j,2*i]=dna["{0}".format(color[j,i])][0]
                dec[j,2*i+1]=dna["{0}".format(color[j,i])][1]
    b_dec=(np.packbits(b_dec,axis=-1))
    g_dec=(np.packbits(g_dec,axis=-1))
    r_dec=(np.packbits(r_dec,axis=-1))
    return b_dec,g_dec,r_dec

def xor_operation_for_decryption(b,g,r,mk):
    m,n = b.shape
    bx=np.chararray((m,n))
    gx=np.chararray((m,n))
    rx=np.chararray((m,n))
    b=b.astype(str)
    g=g.astype(str)
    r=r.astype(str)
    for i in range(0,m):
        for j in range (0,n):
            bx[i,j] = dna["{0}{1}".format(b[i,j],mk[i,j])]
            gx[i,j] = dna["{0}{1}".format(g[i,j],mk[i,j])]
            rx[i,j] = dna["{0}{1}".format(r[i,j],mk[i,j])]

    bx=bx.astype(str)
    gx=gx.astype(str)
    rx=rx.astype(str)
    return bx,gx,rx

def enc_image(b,g,r,iname):
    enc_image = cv2.imread(iname)
    enc_image[:,:,2] = r
    enc_image[:,:,1] = g
    enc_image[:,:,0] = b
    cv2.imwrite(("enc.jpg"), enc_image)
    print("saved ecrypted image as enc.jpg")
    return enc_image

def decrypt(enc_img, fx, fy, fz, fp, Mk, bt, gt, rt):
    r, g, b = split_into_rgb_channels(enc_img)
    p, q = rt.shape
    benc, genc, renc = dna_encode(b, g, r)
    bs, gs, rs = decryption_scrambling(fx, fy, fz, benc, genc, renc)
    bx, rx, gx = xor_operation_for_decryption(bs, gs, rs, Mk)
    blue, green, red = dna_decode(bx, gx, rx)
    green, red = red, green
    recovered_img = np.zeros((p, q, 3), dtype=np.uint8)
    recovered_img[:, :, 0] = red
    recovered_img[:, :, 1] = green
    recovered_img[:, :, 2] = blue

    cv2.imwrite("Recovered.jpg", recovered_img)

#program exec9
if (__name__ == "__main__"):
    file_path = selected_image_path
    print(file_path)
    key,m,n = securekey(file_path)
    update_baker_2d(key,x,y,z)
    blue,green,red=decompose_matrix(file_path)
    blue_e,green_e,red_e=dna_encode(blue,green,red)
    Mk_e = key_matrix_encode(key,blue)
    blue_final, green_final, red_final = xor_operation(blue_e,green_e,red_e,Mk_e)
    x,y,z=gen_chaos_seq(m,n)
    print(x,y,z)
    fx,fy,fz=parallel_sequence_indexing(x,y,z)
    blue_scrambled,green_scrambled,red_scrambled = scramble(fx,fy,fz,blue_final,red_final,green_final)
    b,g,r=dna_decode(blue_scrambled,green_scrambled,red_scrambled)
    enc_img=enc_image(b,g,r,file_path)


    print("decrypting...")
    # decrypt(img,fx,fy,fz,file_path,Mk_e,blue,green,red)
    decrypt(enc_img, fx, fy, fz, file_path, Mk_e, b, g, r)




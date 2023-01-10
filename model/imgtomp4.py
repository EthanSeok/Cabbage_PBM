import re
import os
import numpy as np
import tqdm

path = 'D:/images/cnn_data_2022/paprika/2022.5.28/'
paths = [os.path.join(path , i ) for i in os.listdir(path) if re.search(".JPG$", i )]
## 정렬 작업
imglist = []
for i in paths:
    imglist.append(i)
paths = list(np.sort(imglist))
print(paths)
#len('ims/2/a/2a.2710.png')

pathIn= 'D:/images/cnn_data_2022/paprika/2022.5.28/'
pathOut = './output/test.avi'
fps = 30
import cv2
frame_array = []

for idx , path in tqdm.tqdm(enumerate(paths)):
    img = cv2.imread(path)
    # print(path)
    # print(img)
    height, width, layers = img.shape
    size = (width,height)
    frame_array.append(img)
out = cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
for i in tqdm.tqdm(range(len(frame_array))):
    # writing to a image array
    out.write(frame_array[i])
out.release()
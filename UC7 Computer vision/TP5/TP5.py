import numpy as np
import cv2 as cv
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as Rot
import argparse


def get_info(path_1, path_2, F, Cx, Cy, K1, K2, K3, GetM, GetI):
    print('The image_1 is located at ', path_1)
    print('The image_2 is located at ', path_2)
    print('The F is ', F)
    print('The Cx is ', Cx)
    print('The Cy is ', Cy)
    print('The K1 is ', K1)
    print('The K2 is ', K2)
    print('The K3 is ', K3)
    return(path_1, path_2, float(F), float(Cx), float(Cy), float(K1), float(K2), float(K3), int(GetM), int(GetI))

def preprosImage(path_1, path_2, K, dist_coef):
    img_1 = cv.imread(path_1)
    img_2 = cv.imread(path_2)
    img_1_undis = cv.undistort(img_1, K, dist_coef)
    img_2_undis = cv.undistort(img_2, K, dist_coef)
    gray_1 = cv.cvtColor(img_1_undis,cv.COLOR_BGR2GRAY)
    gray_2 = cv.cvtColor(img_2_undis,cv.COLOR_BGR2GRAY)
    return (gray_1, gray_2)




def getF(gray_1, gray_2):
    sift = cv.xfeatures2d.SIFT_create()
    kp_1, des_1 = sift.detectAndCompute(gray_1,None)
    kp_2, des_2 = sift.detectAndCompute(gray_2,None)
    # create BFMatcher object
    bf = cv.BFMatcher()
    # Match descriptors.
    matches = bf.knnMatch(des_1,des_2, k=2)
    good = []
    pts1 = []
    pts2 = []
    # ratio test as per Lowe's paper
    for i,(m,n) in enumerate(matches):
		    if m.distance < 0.5*n.distance:
			      good.append(m)
			      pts2.append(kp_2[m.trainIdx].pt)
			      pts1.append(kp_1[m.queryIdx].pt)

    pts1 = np.int32(pts1) 
    pts2 = np.int32(pts2)
    F, mask = cv.findFundamentalMat(pts1,pts2,cv.FM_LMEDS)

    # We select only inlier points
    pts1 = pts1[mask.ravel()==1]
    pts2 = pts2[mask.ravel()==1]
    return (F, pts1, pts2)


def drawlines(img1,img2,lines,pts1,pts2):
		''' img1 - image on which we draw the epilines for the points in img2
			  lines - corresponding epilines '''
		r,c = img1.shape
		for r,pt1,pt2 in zip(lines,pts1,pts2):
			  color = tuple(np.random.randint(0,255,3).tolist())
			  x0,y0 = map(int, [0, -r[2]/r[1] ])
			  x1,y1 = map(int, [c, -(r[2]+r[0]*c)/r[1] ])
			  img1 = cv.line(img1, (x0,y0), (x1,y1), color,1)
			  img1 = cv.circle(img1,tuple(pt1),5,color,-1)
			  img2 = cv.circle(img2,tuple(pt2),5,color,-1)
		return img1,img2

def visualize(gray_1 ,gray_2, pts1, pts2):
    # Find epilines corresponding to points in right image (second image) and
    # drawing its lines on left image
    lines1 = cv.computeCorrespondEpilines(pts2.reshape(-1,1,2), 2,F)
    lines1 = lines1.reshape(-1,3)
    img5,img6 = drawlines(gray_1,gray_2,lines1,pts1,pts2)
    # Find epilines corresponding to points in left image (first image) and
		# drawing its lines on right image
    lines2 = cv.computeCorrespondEpilines(pts1.reshape(-1,1,2), 1,F)
    lines2 = lines2.reshape(-1,3)
    img3,img4 = drawlines(gray_2,gray_1,lines2,pts2,pts1)
    plt.subplot(121)
    plt.imshow(img5)
    plt.subplot(122)
    plt.imshow(img3)
    plt.show()


parser = argparse.ArgumentParser(description='Camera motion')
parser.add_argument('--path_1', '-p1', help='The path of the first image', default='data/simulated-pair-of-images/frame0000.jpg')
parser.add_argument('--path_2', '-p2', help='The path of the second image', default='data/simulated-pair-of-images/frame0001.jpg')
parser.add_argument('--F', '-F', help='Focal length', default = 554.254691191187)
parser.add_argument('--Cx', '-Cx', help='Princial point X', default = 320.5)
parser.add_argument('--Cy', '-Cy', help='Principal point Y', default = 240.5)
parser.add_argument('--K1', '-K1', help='The coefficient of order 2', default = 0)
parser.add_argument('--K2', '-K2', help='The coefficient of order 4', default = 0)
parser.add_argument('--K3', '-K3', help='The coefficient of order 6', default = 0)
parser.add_argument('--GetM', '-GM', help='Get the corresponding matrix that you want, 1: get F; 2: get F & E; 3: get F & E & C2Cprime; 4: get F & E & C2Cprime & Euler_angles; 0: nothing', default = 4)
parser.add_argument('--GetI', '-GI', help='Get the two view geometry, 1:get; 2:no', default = 0)
args = parser.parse_args()

if __name__ == '__main__':
    try:
        (path_1, path_2, F, Cx, Cy, K1, K2, K3, GetM, GetI) = get_info(args.path_1, args.path_2, args.F, args.Cx, args.Cy, args.K1, args.K2, args.K3, args.GetM, args.GetI)
    except Exception as e:
        print(e)
    
    #K = np.array([[752.253, 0, 626.696],[0, 752.253, 364.502],[0, 0, 1]])
    #dist_coef = np.array([ 0.0884906, -0.225692, 0, 0, 0.0988774])
    K = np.array([[F, 0, Cx],[0, F, Cy],[0, 0, 1]])
    dist_coef = np.array([K1, K2, 0, 0, K3])
    (gray_1, gray_2) = preprosImage(path_1, path_2, K, dist_coef)
    (F, pts1, pts2)= getF(gray_1, gray_2)
    E=np.dot(np.dot(K.T, F), K)
    pts1 = np.float32(pts1)
    pts2 = np.float32(pts2)
    _, R, t, _ =cv.recoverPose(E,pts1,pts2,K)
    C2Cprime=np.concatenate((R.T,-np.dot(R.T,t)),axis=1)
    r=Rot.from_matrix(C2Cprime[:3,:3])
    Euler_angles=r.as_euler('xyz',degrees=True)
    if(GetM==1):
        print('Matrix F:')
        print(F)
    if(GetM==2):
        print('Matrix F:')
        print(F)
        print('Matrix E:')
        print(E)
    if(GetM==3):
        print('Matrix F:')
        print(F)
        print('Matrix E:')
        print(E)
        print('C2Cprime:')
        print(C2Cprime)
    if(GetM==4):
        print('Matrix F:')
        print(F)
        print('Matrix E:')
        print(E)
        print('C2Cprime:')
        print(C2Cprime)
        print('Euler_angles:')
        print(Euler_angles)
    if(GetI==1):
        visualize(gray_1 ,gray_2, pts1, pts2)







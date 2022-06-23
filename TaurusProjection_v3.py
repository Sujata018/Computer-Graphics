from tkinter import Tk, Canvas, PhotoImage, mainloop, Button, LEFT, RIGHT, Scale, HORIZONTAL, Text
from tkinter.ttk import Frame

from math import sin, cos, radians,sqrt,atan,acos
import numpy as np
from scanLineFill import fillFace

angle_range=370        # range of rotation, in degrees, for drawing circle and rotating it to draw Taurus
angle_step=5          # step in degree, for rotation

viewpoint_x=0          # Initial view point (x,y,z world coordinates, with center of the Taurus at (0,0,0))
viewpoint_y=0
viewpoint_z=40

theta=0                # sperical coordinates for view point
phi=0
r=50

theta_lightSource=0    # sperical coordinates for light source with respect to the world coordinates
phi_lightSource=0
r_lightSource=100

x_lightSource_eye=0    # coordinates of light source with respect to eye coordinates
y_lightSource_eye=0
z_lightSource_eye=0

screen_dist=50         # distance of screen from eye

C=None                 # Canvas to draw the prespective view of the Taurus
img=None               # Image to draw in the canvas 
vertices=[]            # 3-D vertices of the Taurus
vertices_3D_eye=[]
lines=None
faces=[]               # 3-D faces of the Taurus
zColorFaces=[]

inset_size = 100
slide_size = 50
canvas_size= 800

scale=7                # scaling between the inset and main canvas

wireFrame = True
inset=None

iar=1                  # Ambient intensity
iag=1
iab=0
ka=0.3

ipr=1                  # Light source intensity
ipg=1
ipb=0

kd=0.5                 # Diffusion coefficient
ks=0.9                 # Specular coefficient
n=500


def circle(radius,center):
    '''
    Draws a circle with given radius and center, and returns its points
    '''
    global angle_range,angle_step
    a,b,c=center
    for angle in range(0,angle_range,angle_step):
        angle_radians=radians(angle)
        x = a + radius*(cos(angle_radians))
        z = c + radius*(sin(angle_radians))
        vertices.append([x,b,z])
    return vertices

def rotation_matrix(axis, theta):
    '''
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    '''
    axis = np.asarray(axis)
    axis = axis / sqrt(np.dot(axis, axis))
    a = cos(theta / 2.0)
    b, c, d = -axis * sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def drawTaurus(axis,circle_vertices):
    '''
    Rotate circle_vertices around axis to form the points for a taurus
    '''
    global angle_range,angle_step

    circle_vertices_rotated = []
    for angle in range(0,angle_range,angle_step):
        angle_radians=radians(angle)
        for vertex in circle_vertices:
            vertex_rotated=np.dot(rotation_matrix(axis,angle_radians),vertex)
            circle_vertices_rotated.append(vertex_rotated)
    return circle_vertices_rotated

def drawFaces(num_points_row,num_cols):
    '''
    Determine the faces of the taurus
    '''
    Tfaces =[]
    for col in range(num_cols-1):
        for i in range(num_points_row-1):  
            pt1=i+col*num_points_row
            pt2=i+(col+1)*num_points_row
            face=[pt1,pt2,pt2+1,pt1+1]  
            Tfaces.append(face)
    return Tfaces

def drawPatchModel():
    '''
    color patches
    '''
    global wireFrame,img

    wireFrame=False
    img=refresh_canvas(img)

def drawWireFrameModel():
    '''
    draw wireframes
    '''
    global wireFrame,img

    wireFrame=True
    img=refresh_canvas(img)
    return img

def determineZofFace(face):
    global vertices_3D_eye
    z1=vertices_3D_eye[face[0]][2]
    z2=vertices_3D_eye[face[1]][2]
    z3=vertices_3D_eye[face[2]][2]
    z4=vertices_3D_eye[face[3]][2]
    return (z1+z2+z3+z4)/4

def drawWireFrame():
    '''
    Draws the outlines of the faces to create a wireframe model, in the main canvas
    '''
    global lines,faces,C,scale,vertices_2D,img
    
    shift=canvas_size//2                                 # shifting origin to mid point of canvas for full view of taurus in positive coordinates

    zfaces=[]
    for face in faces:
        z=determineZofFace(face)
        zfaces.append((face[0],face[1],face[2],face[3],z))

    zfaces.sort(key=lambda x:x[4],reverse=True)
    for face in zfaces:
        x1,y1=vertices_2D[face[0]]*scale+shift
        x2,y2=vertices_2D[face[1]]*scale+shift
        x3,y3=vertices_2D[face[2]]*scale+shift
        x4,y4=vertices_2D[face[3]]*scale+shift

        C.create_polygon(x1,y1,x2,y2,x3,y3,x4,y4,fill="#FFFFFF",outline="#000000")
        
    return img

def normalise(x,y,z):
    
    magn=sqrt(x**2+y**2+z**2)
    return x/magn,y/magn,z/magn

def pIllum(face,patchColor):

    global iar,iag,iab,ka,ipr,ipg,ipb,kd,ks,x_lightSource_eye,y_lightSource_eye,z_lightSource_eye,n

    ## Extract RGB components of patch color
    pr=int(patchColor[1:3],16)/256
    pb=int(patchColor[5:],16)/256
    pg=int(patchColor[3:5],16)/256
    
    ## Add ambient component
    pr += ka*iar
    pg += ka*iag
    pb += ka*iab

    ## Get middle point of the patch
    x1,y1,z1=vertices_3D_eye[face[0]]
    x2,y2,z2=vertices_3D_eye[face[1]]
    x3,y3,z3=vertices_3D_eye[face[2]]
    ix=(x1+x3)//2
    iy=(y1+y3)//2
    iz=(z1+z3)//2
    ix,iy,iz=normalise(ix,iy,iz)
    
    ## Determine L vector
    lx=x_lightSource_eye-ix
    ly=y_lightSource_eye-iy
    lz=z_lightSource_eye-iz
    lx,ly,lz=normalise(lx,ly,lz)

    ## Determine attenuation factor
    fatt=1/(lx**2+ly**2+lz**2)
    
    ## Determine A, B, C of the plane for plane equation Ax+By+Cz+D=0
    A=y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2)
    B=x1*(z2-z3)+x2*(z3-z1)+x3*(z1-z2)
    C=x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
    A,B,C=normalise(A,B,C)

    ## Determine cosTheta  (N . L)
    cosTheta=lx*A+ly*B+lz*C

    ## Add diffusion component
    if cosTheta > 0:
        pr += fatt*ipr*kd*cosTheta
        pg += fatt*ipg*kd*cosTheta
        pb += fatt*ipb*kd*cosTheta

    ## Determine R vector
    rx=2*A*cosTheta-lx
    ry=2*B*cosTheta-ly
    rz=2*C*cosTheta-lz
    rx,ry,rz=normalise(rx,ry,rz)

    ## Determine cosAlpha  (V . R)
    cosAlpha=(-1)*(rx*ix+ry*iy+rz*iz)

    ## Add spacular component
    if cosAlpha > 0:
        pr += fatt*ipr*ks*(cosAlpha**n)
        pg += fatt*ipg*ks*(cosAlpha**n)
        pb += fatt*ipb*ks*(cosAlpha**n)

    ## Convert RGB to hex format for Tkinter
    if pr > 1:
        pr=1
    if pg > 1:
        pg=1
    if pb > 1:
        pb=1
    
    pr=int(pr*255)
    pg=int(pg*255)
    pb=int(pb*255)
    
    hexa_r=hex(int(pr))[2:].zfill(2)
    hexa_g=hex(int(pg))[2:].zfill(2)
    hexa_b=hex(int(pb))[2:].zfill(2)
    illumColor='#'+ hexa_r + hexa_g + hexa_b
        
    return illumColor

def drawPatches():
    '''
    Fills the faces to create a patch model, in the main canvas
    '''
    global faces,C,scale,vertices_2D,img,zColorFaces,theta,phi,r,x_lightSource_eye,y_lightSource_eye,z_lightSource_eye,r_lightSource,theta_lightSource,phi_lightSource
    
    shift=canvas_size//2      # shifting origin to mid point of canvas for full view of taurus in positive coordinates

    # Convert light source world coordinates to eye coordinates
    x=r_lightSource*sin(phi_lightSource)*cos(theta_lightSource)
    y=r_lightSource*sin(phi_lightSource)*sin(theta_lightSource)
    z=r_lightSource*cos(phi_lightSource)

    lightSource_eye=get3DeyeFrom3DWorld([[x,y,z]])
    x_lightSource_eye,y_lightSource_eye,z_lightSource_eye=lightSource_eye[0]

    patchColor="#00FF00"
    dir = 1
    zfaces=[]
    for face in faces:
        z=determineZofFace(face)
        zfaces.append((face[0],face[1],face[2],face[3],z))

    zfaces.sort(key=lambda x:x[4],reverse=True)
    for face in zfaces:
        illumColor=pIllum(face,patchColor)
        x1,y1=vertices_2D[face[0]]*scale+shift
        x2,y2=vertices_2D[face[1]]*scale+shift
        x3,y3=vertices_2D[face[2]]*scale+shift
        x4,y4=vertices_2D[face[3]]*scale+shift
        C.create_polygon(x1,y1,x2,y2,x3,y3,x4,y4,fill=illumColor)
    return img

def get3DeyeFrom3DWorld(vertices):
    '''
    Convert from world coordinates to view point coordinates
    '''
    global theta,phi,r

    V=np.array([(-1)*sin(theta),(-1)*cos(phi)*cos(theta),(-1)*sin(phi)*cos(theta),0,cos(theta),(-1)*cos(phi)*sin(theta),(-1)*sin(phi)*sin(theta),0,0,sin(phi),(-1)*cos(phi),0,0,0,r,1 ]).reshape([4,4])

    vertices_homo=np.array(vertices)                                            # calculating projected points in 2D using homogeneous coordinates through matrix multiplication
    vertices_homo_3D_world=np.hstack([vertices_homo,np.ones(len(vertices)).reshape([len(vertices),1])])

    vertices_eye=np.dot(vertices_homo_3D_world,V)
    return vertices_eye[:,:3]

def get2DperspectiveFrom3DEye():
    '''
    Get 2D perspective coordinates from 3D eye coordinates
    '''
    global screen_dist, vertices_2D, vertices_3D_eye

    P=np.array([screen_dist,0,0,0,screen_dist,0,0,0,1]).reshape([3,3])
    vertices_perspective = np.dot(vertices_3D_eye,P)

    vertices_2D=np.zeros(vertices_perspective.shape)
    vertices_2D=vertices_2D[:,:2]

    vertices_2D[:,0]=vertices_perspective[:,0]/vertices_perspective[:,2]
    vertices_2D[:,1]=vertices_perspective[:,1]/vertices_perspective[:,2]

    return

def refresh_canvas(img):
    '''
    Draw projection of Taurus in main canvas
    '''
    global r,theta,phi,vertices,wireFrame,vertices_3D_eye

    C.delete('all')     # Clear the canvas
    img=PhotoImage(height=canvas_size,width=canvas_size)
    C.create_image((canvas_size//2,canvas_size//2),image=img, state="normal")

    vertices_3D_eye=get3DeyeFrom3DWorld(vertices)

    get2DperspectiveFrom3DEye()
    
    if wireFrame == True:
        img=drawWireFrame()
    else:
        img=drawPatches()
    return img

def change_theta_phi_r(action):
    '''
    Used in button function to receive sperical coordinates of the view point
    '''
    global theta,phi,r,img,inset
    
    theta=sli_theta.get()
    phi=sli_phi.get()
    r=sli_r.get()

    img=refresh_canvas(img)
    drawInset()
    
    return img

def change_theta_phi_r_light_source(action):
    '''
    Used in button function to receive sperical coordinates of the light source
    '''
    global r_lightSource,theta_lightSource,phi_lightSource,img,inset
    
    theta_lightSource=sli_lightSource_theta.get()
    phi_lightSource=sli_lightSource_phi.get()
    r_lightSource=sli_lightSource_r.get()

    img=refresh_canvas(img)
    drawInset()
    
    return img

def drawInset():

    global inset,theta,phi,r,theta_lightSource,phi_lightSource,r_lightSource

    inset.delete('all')     # Clear the canvas
    
    inset.create_text(80,70,text='r=200')
    inset.create_oval(25,25,75,75,width=1)
    inset.create_line(25,50,75,50,width=1)
    inset.create_oval(30,45,70,55,fill="blue")
    inset.create_oval(33,48,67,52,fill="white")

    inset.create_oval(65,3,71,9,fill="red")
    inset.create_text(85,6,text='eye')
    inset.create_oval(65,15,71,21,fill="yellow")
    inset.create_text(85,18,text='light')

    x=r*cos(theta)*sin(phi)//8+50
    z=50-r*cos(phi)//8

    inset.create_oval(x-3,z-3,x+3,z+3,fill="red")

    x_ls=r_lightSource*cos(theta_lightSource)*sin(phi_lightSource)//8+50
    z_ls=50 - r_lightSource*cos(phi_lightSource)//8

    inset.create_oval(x_ls-3,z_ls-3,x_ls+3,z_ls+3,fill="yellow")
    
    return

if __name__=='__main__':

    # Create points for the circle
    radius=10
    a=20         # x,y,z components of the center
    b=0
    c=0

    num_points_row=int(angle_range/angle_step)   # Total number of points in a row

    circle_vertices=circle(radius,(a,b,c))       # Get points for a circle

    axis = [0, 0, 1]                             # z-axis

    num_cols=int(angle_range//angle_step)        # number of circles in a taurus

    vertices=drawTaurus(axis,circle_vertices)    # Get all points for the taurus

    ## Create faces of the Taurus

    faces=drawFaces(num_points_row,num_cols)

    ## Show Taurus using Tkinter canvases

    top=Tk()                                                         # Create an instance of tkinter frame
    inset=Canvas(top,bg="white",height=inset_size,width=inset_size)  # Define inset canvas in white
    
    sli_theta=Scale(top,label="Eye Theta",from_=0, to =359, orient=HORIZONTAL)
    sli_phi=Scale(top,label="Eye Phi",from_=0, to =90, orient=HORIZONTAL)
    sli_r=Scale(top,label="Eye r",from_=0, to =200, orient=HORIZONTAL)

    sli_theta.set(theta)
    sli_phi.set(phi)
    sli_r.set(r)
    
    sli_lightSource_theta=Scale(top,label="Light source Theta",from_=0, to =359, orient=HORIZONTAL)
    sli_lightSource_phi=Scale(top,label="Light source Phi",from_=0, to =90, orient=HORIZONTAL)
    sli_lightSource_r=Scale(top,label="Light source r",from_=0, to =200, orient=HORIZONTAL)

    sli_lightSource_theta.set(theta_lightSource)
    sli_lightSource_phi.set(phi_lightSource)
    sli_lightSource_r.set(r_lightSource)
    
    button_wireFrame = Button(top, text="WireFrame", command=drawWireFrameModel)
    button_Patch = Button(top, text="Patches", command=drawPatchModel)
    drawInset()    

    C=Canvas(top,bg="white",height=canvas_size,width=canvas_size)    # Define main canvas C, in white

    C.pack()

    inset.place(x=100,y=5)
    sli_theta.place(x=250,y=1)
    sli_phi.place(x=400,y=1)
    sli_r.place(x=550,y=1)
 
    sli_lightSource_theta.place(x=250,y=50)
    sli_lightSource_phi.place(x=400,y=50)
    sli_lightSource_r.place(x=550,y=50)
    button_wireFrame.place(x=5,y=5)
    button_Patch.place(x=5,y=55)
    C.place(x=5,y=101)

    img=PhotoImage(height=canvas_size,width=canvas_size)
    C.create_image((canvas_size//2,canvas_size//2),image=img, state="normal")

    img=refresh_canvas(img)                                          # Draw initial taurus from viewpoint on z-axis

    sli_theta.bind('<ButtonRelease>',change_theta_phi_r) 
    sli_phi.bind('<ButtonRelease>',change_theta_phi_r)    
    sli_r.bind('<ButtonRelease>',change_theta_phi_r)         

    sli_lightSource_theta.bind('<ButtonRelease>',change_theta_phi_r_light_source) 
    sli_lightSource_phi.bind('<ButtonRelease>',change_theta_phi_r_light_source)    
    sli_lightSource_r.bind('<ButtonRelease>',change_theta_phi_r_light_source)         
                                                             
    top.mainloop()

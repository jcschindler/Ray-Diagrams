import math
import numpy as np
import matplotlib.pyplot as plt

INDEX_AIR = 1.000

def mag(a):
    return math.sqrt(np.sum(a*a))

def intersection(a,b):  # a-light pos vec, b - element pos vec (lens ->x[1,:,:])
    c = a[1]-a[0] #define direction vecs
    d = b[1]-b[0]
    w = a[0]- b[0]
    s = ((d[1]*w[0]) - (d[0]*w[1]))/((d[0]*c[1])-(d[1]*c[0]))
    return s*c + a[0]


class light:
    def __init__(self,lens,distance,number_of_rays):
        self.loc = [-distance,lens.height/2]
        self.nRays = number_of_rays
        self.coresponding_el = np.zeros(self.nRays)

    def initialize_rays(self,lens): # make light at [0,0] to reduce calc time in future
        y = np.linspace(lens.height/(self.nRays+2),(self.nRays+1)*lens.height/(self.nRays+2),self.nRays)
        eRays =  np.zeros([3,self.nRays,2,2])#ray segments [preLens-Lens-postLens,nrays,1-2,x-y]
        eRays[0,:,0,:] = self.loc
        for i in range(self.nRays):
            eRays[0,i,1,:] = [0,y[i]]
        return eRays

    def first_contact(self,lens,x,eRays):
        y = np.linspace(0,lens.height,lens.N)
        for i in range(self.nRays):
            which_el = -2
            which_el_old = -100
            while(abs(which_el - which_el_old) >= 1):
                which_el_old = which_el
                which_el = max(max(np.nonzero(eRays[0,i,1,1]>y)))  # change 0 for later collision
                # print(which_el)
                eRays[0,i,1,:] = intersection(eRays[0,i],x[which_el])
                self.coresponding_el[i] = which_el
        eRays[1,:,0,:] = eRays[0,:,1,:]
        return eRays

    def ray_from_theta(self,which_time,theta2,eRays,x):

        if which_time == 0:
            for i in range(self.nRays):
                x1 = math.cos(theta2[i]) #component in direction of normal
                y1 = math.sin(theta2[i])   #component in dir of elements
                element_dir = x[int(self.coresponding_el[i]),1,:]-x[int(self.coresponding_el[i]),0,:]
                element_dir = element_dir/mag(element_dir)
                # print(mag(element_dir))
                n_dir = [-element_dir[1],element_dir[0]]
                eRays[1,i,1,:] = [(element_dir[0]*y1)+(n_dir[0]*x1)+eRays[1,i,0,0] , (element_dir[1]*y1)+(n_dir[1]*x1)+eRays[1,i,0,1] ]
                eRays[2,i,0,:] = eRays[1,i,1,:]
        if which_time == 1:
            scale = 50
            for i in range(self.nRays):
                x1 = scale*math.cos(theta2[i])
                y1 = scale*math.sin(theta2[i])
                element_dir = x[int(self.coresponding_el[i]),1,:]-x[int(self.coresponding_el[i]),0,:]
                element_dir = element_dir/mag(element_dir)
                n_dir = [-element_dir[1],element_dir[0]]
                eRays[2,i,1,:] = [-((element_dir[0]*y1)+(n_dir[0]*x1)+eRays[2,i,0,0]),-(element_dir[1]*y1)-(n_dir[1]*x1) +eRays[2,i,0,1]]


    # def second_contact(self,lens,x,eRays):
    #     y = np.linspace(0,lens.height,lens.N)
    #     for i in range(self.nRays):
    #         which_el = -2
    #         which_el_old = -100
    #         while(abs(which_el - which_el_old) >= 1):
    #             which_el_old = which_el
    #             which_el = lens.N-1+max(max(np.nonzero(eRays[1,i,1,1]>y)))  # change 0 for later collision
    #             # print(which_el)
    #             eRays[1,i,1,:] = intersection(eRays[1,i],x[which_el])
    #             self.coresponding_el[i] = which_el
    #     eRays[2,:,0,:] = eRays[1,:,1,:]
    #     return eRays

    def second_contact(self,lens,x,eRays):
        y = np.linspace(0,lens.height,lens.N)
        for i in range(self.nRays):
            which_el = -2
            which_el_old = -100
            while(abs(which_el - which_el_old) >= 1):
                which_el_old = which_el
                # print(np.shape(eRays[1,i,1,1]))
                if max(np.nonzero(eRays[1,i,1,1]>y)).size!=0:
                    which_el = lens.N-1+max(max(np.nonzero(eRays[1,i,1,1]>y)))  # change 0 for later collision
                else:
                    which_el = 2*lens.N-2
                if which_el >= 2*lens.N-2:
                    eRays[1,i,:,:] = np.zeros([2,2])
                    eRays[0,i,:,:] = np.zeros([2,2])
                    eRays[2,i,:,:] = np.zeros([2,2])
                    which_el_old = which_el
                else:
                    eRays[1,i,1,:] = intersection(eRays[1,i],x[which_el])
                    self.coresponding_el[i] = which_el

        eRays[2,:,0,:] = eRays[1,:,1,:]
        return eRays



class lens:
    def __init__(self,index,number_of_elements,height,thickness,symmetry):
        self.index = index
        self.N = number_of_elements + 1 #odd number please (num of points, num el = 2(N-1))
        self.thickness = thickness
        self.height = height
        # self.shape = 1 #half circle
        self.symmetry = symmetry # [-1:forward , 0:symmetric , 1:reverse

    def get_widths(self,shape,a4,a6,a8,k,R_over_H,R):
        if shape == 'circle':
            angles = np.linspace(-math.pi/2,math.pi/2,self.N)
            y = np.linspace(0,self.height,self.N)
            return np.cos(angles)+self.thickness,y
        elif shape == 'circle_divergent':
            angles = np.linspace(-math.pi/2,math.pi/2,self.N)
            y = np.linspace(0,self.height,self.N)
            return 2-np.cos(angles)+self.thickness,y
        elif shape == 'triangle':
            base = self.height/2
            y = np.linspace(0,self.height,self.N)
            return base-(y/2),y
        elif shape == 'sine':
            w = 2*2*math.pi/self.height
            y = np.linspace(0,self.height,self.N)
            return .1*np.sin(w*(y-self.height/2)**2-1.6)+self.thickness,y
        elif shape == 'combo':
            angles = np.linspace(-math.pi/2,math.pi/2,self.N)
            w = 1*2*math.pi/self.height
            y = np.linspace(0,self.height,self.N)
            return .10*np.cos(w*y)+self.thickness+1*np.cos(angles),y
        elif shape == 'asphere':
            R = self.height*R_over_H
            y = np.linspace(0,self.height,self.N)
            r = y - (self.height/2)
            z = (r**2/(R*(1+np.sqrt(1-((1+k)*(r**2/R**2))))))+(a4*r**4)+(a6*r**6)
            return R-z+self.thickness,y
        elif shape == 'diode_laser':
            R = R
            y = np.linspace(0,self.height,self.N)
            r = y - (self.height/2)
            z = (r**2/(R*(1+np.sqrt(1-((1+k)*(r**2/R**2))))))+(a4*r**4)+(a6*r**6)+(a8*r**8)
            return R-z,y

    def assemble_lens(self,lens_shape,a4 = 0,a6=0,a8 = 0,k=-1,R_over_H = 3/4,R = 0):
        widths,y = self.get_widths(shape=lens_shape,a4 = a4,a6=a6,a8 = a8,k=k,R_over_H = R_over_H,R=R)
        # x =
        if self.symmetry == 0:
            # this separates out elements ---v
            x=np.zeros([2*(self.N-1),2,2])
            for i in range((self.N-1)): #start bottom CCW, element definitions
                x[i,0,:] = [-widths[i]/2,y[i]]
                x[i,1,:] = [-widths[i+1]/2,y[i+1]]
                x[i+self.N-1,0,:] = [widths[i]/2,y[i]]
                x[i+self.N-1,1,:] = [widths[i+1]/2,y[i+1]]
            return x #elements for plotting
        elif self.symmetry == -1:
            x=np.zeros([2*(self.N-1),2,2])
            for i in range((self.N-1)): #start bottom CCW, element definitions
                x[i,0,:] = [-widths[i],y[i]]
                x[i,1,:] = [-widths[i+1],y[i+1]]
                x[i+self.N-1,0,:] = [0,y[i]]
                x[i+self.N-1,1,:] = [0,y[i+1]]
            return x #elements for plotting
        elif self.symmetry == 1:
            x=np.zeros([2*(self.N-1),2,2])
            for i in range((self.N-1)): #start bottom CCW, element definitions
                x[i+self.N-1,0,:] = [widths[i],y[i]]
                x[i+self.N-1,1,:] = [widths[i+1],y[i+1]]
                x[i,0,:] = [0,y[i]]
                x[i,1,:] = [0,y[i+1]]
            return x #elements for plotting

    def snell(self,which_time,light,eRays,x,n_other):
        theta1 = np.zeros(light.nRays)
        evRays = eRays[which_time,:,1,:]-eRays[which_time,:,0,:]
        xv = x[:,1,:]-x[:,0,:]
        for i in range(light.nRays):
            theta1[i] = math.acos(np.sum(evRays[i]*xv[int(light.coresponding_el[i])])/(mag(xv[int(light.coresponding_el[i])])*mag(evRays[i])))-math.pi/2
        if which_time==0:
            theta2 = np.arcsin(n_other*np.sin(theta1)/self.index)
            light.ray_from_theta(which_time,theta2,eRays,x)
        elif  which_time == 1:
            theta2 = np.arcsin(self.index*np.sin(theta1)/n_other)
            light.ray_from_theta(which_time,theta2,eRays,x)

def run(lens,light,x_lens,INDEX_AIR):
    e_ray = light.initialize_rays(lens)
    light.first_contact(lens,x_lens,e_ray)
    lens.snell(0,light,e_ray,x_lens,INDEX_AIR)
    light.second_contact(lens,x_lens,e_ray)
    lens.snell(1,light,e_ray,x_lens,INDEX_AIR)
    plt.axes([.1,.1,.8,.8], aspect=1.)
    for i in range(2*lens.N-2):
        plt.plot(x_lens[i,:,0],x_lens[i,:,1],color=[.25,.25,.25],linewidth=.9)
    plt.plot([x_lens[0,0,0],x_lens[lens.N-1,0,0]],[x_lens[0,0,1],x_lens[lens.N-1,0,1]],color=[.25,.25,.25],linewidth=.9)
    plt.plot([x_lens[lens.N-2,1,0],x_lens[-1,1,0]],[x_lens[lens.N-2,1,1],x_lens[-1,1,1]],color=[.25,.25,.25],linewidth=.9)
    for section in range(3):
        for i in range(light.nRays):
            plt.plot(e_ray[section,i,:,0],e_ray[section,i,:,1],'k',linewidth=.5)
    plt.axis([-10,30,-lens.height,lens.height*2])
    plt.show()
#
# asphere
# lens1 = lens(index=1.75,number_of_elements=500,height=10,thickness=-5,symmetry=0)
# lamp = light(lens1,distance=10000000,number_of_rays=15)
# x_lens = lens1.assemble_lens(lens_shape = 'asphere',a4 = -0.0002, a6=-.0000, k=-1, R_over_H = 3/4)
# # x_lens = lens1.assemble_lens(lens_shape = 'asphere',a4 = -0.0002, a6=-.0001, k=-1, R_over_H = 3/4)
# run(lens1,lamp,x_lens,INDEX_AIR)

# spherical
# lens2 = lens(index=1.75,number_of_elements=100,height=10,thickness=0,symmetry=0)
# lamp2 = light(lens2,distance=100000,number_of_rays=15)
# x_lens2 = lens2.assemble_lens(lens_shape = 'circle')
# run(lens2,lamp2,x_lens2,INDEX_AIR)

# spherical divergent
# lens2 = lens(index=1.75,number_of_elements=100,height=10,thickness=0,symmetry=0)
# lamp2 = light(lens2,distance=100000,number_of_rays=15)
# x_lens2 = lens2.assemble_lens(lens_shape = 'circle_divergent')
# run(lens2,lamp2,x_lens2,INDEX_AIR)

# triangle
# lens2 = lens(index=1.75,number_of_elements=1,height=10,thickness=0,symmetry=0)
# lamp2 = light(lens2,distance=20,number_of_rays=10)
# x_lens2 = lens2.assemble_lens(lens_shape = 'triangle')
# run(lens2,lamp2,x_lens2,INDEX_AIR)

# sinusoidal
# lens2 = lens(index=1.75,number_of_elements=500,height=10,thickness=1,symmetry=0)
# lamp2 = light(lens2,distance=10,number_of_rays=101)
# x_lens2 = lens2.assemble_lens(lens_shape = 'sine')
# run(lens2,lamp2,x_lens2,INDEX_AIR)

# cosine + shperical
lens2 = lens(index=1.75,number_of_elements=500,height=10,thickness=1,symmetry=0)
lamp2 = light(lens2,distance=10000,number_of_rays=21)
x_lens2 = lens2.assemble_lens(lens_shape = 'combo')
run(lens2,lamp2,x_lens2,INDEX_AIR)

# ## diode laser system
# lens2 = lens(index=2.604,number_of_elements=1000,height=4,thickness=0,symmetry=1)
# lamp2 = light(lens2,distance=1.75,number_of_rays=2)
# x_lens2 = lens2.assemble_lens(lens_shape = 'diode_laser',a4 = -0.004027423,a6=-0.000850078,a8 = 0.0000045098755,k=-0.06201888,R = 3)
# run(lens2,lamp2,x_lens2,INDEX_AIR)









# e_ray = lamp.initialize_rays(lens1)
# lamp.first_contact(lens1,x_lens,e_ray)
# lens1.snell(0,lamp,e_ray,x_lens,INDEX_AIR)
# lamp.second_contact(lens1,x_lens,e_ray)
# lens1.snell(1,lamp,e_ray,x_lens,INDEX_AIR)
# plt.axes([.1,.1,.8,.8], aspect=1.)
# for i in range(2*lens1.N-2):
#     plt.plot(x_lens[i,:,0],x_lens[i,:,1],color=[.25,.25,.25],linewidth=.9)
# for section in range(3):
#     for i in range(lamp.nRays):
#         plt.plot(e_ray[section,i,:,0],e_ray[section,i,:,1],'k',linewidth=.5)
# plt.axis([-10,30,-lens1.height,lens1.height*2])
# plt.show()

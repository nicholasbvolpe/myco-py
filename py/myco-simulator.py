import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
import time
from math import sin, cos
from random import uniform


def normalize(vector):
    '''
    Normalizes the given vector, used to make direction vectors
    for the field
    '''
    return vector/norm(vector)


class section():
    '''
    A structure defining the attributes of a section of mycelium
    '''

    def __init__(self, start_point, end_point, growth_vec=[0, 0, 0]):
        # Start point of section (3 dimensional point)
        self.f = np.array(start_point)
        # End point of section (3 dimensional point)
        self.t = np.array(end_point)
        self.age = 0  # Age is incremented every timestep. Used to prevent excessive branches
        self.g = np.array(growth_vec)
        self.branches = list()
        # A list of points that allow a section to 'curve' from tropisms
        # Not currently implmeneted
        # self.sub_segments = list()


class MyceliumSimulator():
    def __init__(self):
        self.length = 0  # Total length of mycelium
        self.n = 0  # Number of hyphal sections
        self.sections = []  # Stores all sections

    def add_sec(self, s):
        '''
        Add a mycelial section to the simulator
        '''
        self.sections.append(s)
        self.n += 1

    def gen_plot(self, param_str):
        '''
        Plot all sections as well as their
        '''
        if self.sections == []:
            print("There are no sections to plot")
        else:
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            for section in self.sections:
                # Plot sections as well as their starting coords, and plots end coords
                # of the sections without branches (tips)
                if section.branches == []:
                    ax.scatter(section.t[0], section.t[1],
                               section.t[2], color='red', s=3)
                ax.scatter(section.f[0], section.f[1],
                           section.f[2], color='white', s=3)
                ax.plot([section.f[0], section.t[0]], [section.f[1], section.t[1]], [
                        section.f[2], section.t[2]], color="brown")
            if param_str != "":
                plt.title(param_str)
            plt.show()

    def hyphal_avoid_field(self, p, index, secLen):
        N_Sp = 0  # Stores the total field
        D = []  # Stores vectors to p from half-way points of sections in L
        m = 100
        for x in range(secLen):
            if index != x:
                P_m = []
                # Stores field values for m points
                fields = 0 
                # Creating vector of m points
                P_m.extend(zip(np.linspace(self.sections[x].f[0], self.sections[x].t[0], m), np.linspace(
                    self.sections[x].f[1], self.sections[x].t[1], m), np.linspace(self.sections[x].f[2], self.sections[x].t[2], m)))
                for k in P_m:
                    # Finding field value for each of the m points
                    val = (norm(p-np.array(k))**2)
                    fields += (1/val)
                # The field for this section is added to N_Sp
                N_Sp += (fields/m)
                #d_point = [(self.sections[x].f[0]+self.sections[x].t[0])*0.5, (self.sections[x].f[1]+self.sections[x].t[1])
                #        * 0.5, (self.sections[x].f[2]+self.sections[x].t[2])*0.5]
                #print(d_point)
                D.append(p-np.array(self.sections[x].t))
        # The total field at point p is equal to the sum of the fields on p from every section
        # The direction of the field is the sum of all vectors in D
        dlen = len(D)
        #if dlen == 0:
        #    return np.nan, np.nan
        d = np.sum(D, axis=0)
        unit_vec = d/norm(d)
        return N_Sp, unit_vec

    def neighbors(self, j, dist):
        '''
        Counts how many sections are within a distance equal to dist from section j.
        '''
        neighbor_count = 0
        for i in range(len(self.sections)):
            if j != i and not np.array_equal(self.sections[i].t, self.sections[j].f):
                if norm(self.sections[j].t-self.sections[i].t) < dist or norm(self.sections[j].t-self.sections[i].f) < dist:
                    neighbor_count += 1
        return neighbor_count

    def run_neighbor(self, ITERATIONS, max_dist=15, max_branch_dist=15, max_nbrs=15, max_branch_nbrs=3, branch_probability=0.40, a=1):
        first = time.perf_counter()
        if self.n <= 0:
            print("No sections")
        for step_i in range(ITERATIONS):
            for j in range(len(self.sections)):
                if len(self.sections[j].branches) == 0:
                    if self.neighbors(j, max_dist) < max_nbrs:
                        # Growth step
                        step_length = a * self.sections[j].g
                        self.length += a
                        self.sections[j].t = self.sections[j].t + step_length

                        # Growth vector is not updated in the simple neighbors mode of simulation

                        # Branching Step
                            # Branch if there are less than a certain number of sections nearby
                            # and if branch_probability is satisfied
                            # Age is restricted so that branches aren't created one after another
                        branchDecision = uniform(0, 1)
                        cur_age = self.sections[j].age
                        if branchDecision < branch_probability and self.neighbors(j, max_branch_dist) < max_branch_nbrs:
                            # Creating random angles for branching
                            rot1 = [uniform(-60, 61), uniform(-60, 61)]
                            rot2 = [uniform(-60, 61), uniform(-60, 61)]
                            # Use random angles to create a scipy rotation object
                            r1 = R.from_euler('ZY', rot1, degrees=True)
                            r2 = R.from_euler('ZY', rot2, degrees=True)
                            # The rotation is applied to the section's growth vector
                            new_g_vec1 = r1.apply(self.sections[j].g)
                            new_g_vec2 = r1.apply(self.sections[j].g)
                            # Create branches and add them to section j's branch list
                            # One branch maintains the parent section's growth vector to reflect
                            # the biology of mycelium
                            branch1 = section(self.sections[j].t, self.sections[j].t, new_g_vec1)
                            branch2 = section(self.sections[j].t, self.sections[j].t, new_g_vec2)
                            self.sections[j].branches.append(branch1)
                            self.sections[j].branches.append(branch2)
                            # Adds branches to the simulator's list of all sections
                            self.add_sec(branch1)
                            self.add_sec(branch2)
                self.sections[j].age += 1
        second = time.perf_counter()
        print("Time to run simulation: ", second-first)
        print("Number of sections: ", self.n)
        print("Total length of mycelium: ", self.length)
        param_str = "Iterations: "+str(ITERATIONS)+" Branch_prob: "+str(branch_probability)+" Neighbor limit: "+str(
            max_nbrs)+" Neighbor dist: "+str(max_dist)+"\nBranch Neighbor limit: "+str(max_branch_nbrs)+" Branch Neighbor dist: "+str(max_branch_dist)
        sim.gen_plot(param_str)

    def run_tropisms(self, ITERATIONS, limit_f1=0.1, branch_probability=0.40, a=1, k=0.5):
        first = time.perf_counter()
        if self.n <= 0: print("No sections")
        for step_i in range(ITERATIONS):
            secLen = len(self.sections)
            for j in range(secLen):
                # Hyphal avoidance field
                # Sections that have branched are unable to grow or branch further
                if len(self.sections[j].branches) == 0:
                    # Growth step
                    step_length = a * self.sections[j].g
                    self.length += a
                    self.sections[j].t = self.sections[j].t + step_length

                    # Each field returns scalar value of field and unit vector
                    # of field on current section's end point (self.sections[j].t)
                    f1, f1u = self.hyphal_avoid_field(self.sections[j].t, j, secLen)

                    vs = f1*f1u
                    vs_n = normalize(vs)

                    # Calculate new growth vector
                    new_growth = k*self.sections[j].g+(1-k)*vs_n
                    self.sections[j].g = normalize(new_growth)

                    # Branching Step
                    cur_age = self.sections[j].age 
                    if cur_age > 0 and uniform(0, 1) < branch_probability and f1 < limit_f1:
                        # Creating random growth directions for second branch
                        rot2 = [uniform(-60, 61), uniform(-60, 61)]
                        # Use random angles to create a scipy rotation object
                        r2 = R.from_euler('ZY', rot2, degrees=True)
                        new_g_vec2 = r2.apply(self.sections[j].g)
                        # Create and add branches to section j's branch list
                        branch1 = section(self.sections[j].t, self.sections[j].t, new_growth)
                        branch2 = section(self.sections[j].t, self.sections[j].t, new_g_vec2)
                        self.sections[j].branches.append(branch1)
                        self.sections[j].branches.append(branch2)
                        self.add_sec(branch1)
                        self.add_sec(branch2)

                self.sections[j].age += 1
        second = time.perf_counter()
        print("Time to run simulation: ", second-first)
        print("Number of sections: ", self.n)
        print("Total length of mycelium: ", self.length)
        param_str = "Iterations: "+str(ITERATIONS)+" Branch_prob: "+str(
            branch_probability)+" Growth Rate: "+str(a)+" Persistence Factor: "+str(k)
        sim.gen_plot(param_str)
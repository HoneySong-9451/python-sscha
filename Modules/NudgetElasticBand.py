from __future__ import print_function
from __future__ import division

import cellconstructor as CC 
import cellconstructor.Phonons 

import sys, os

import numpy as np 
import Cluster, SchaMinimizer, Ensemble

import multiprocessing 

class NEB:
    """
    NUDGET ELASTIC BAND
    ===================

    This class is done for identifying the transition path between two states.
    It works by creating many replica of the system and attaching springs
    in between them, keeping fixed the starting and final position.
    The optimiziation is performed by nudging, i.e. the perpendicular direction of the
    band will be optimized according to the potential, while the longitudinal direction is
    optimized according to the spring (this avoids cut of cross errors).
    """
    

    def __init__(self):
        self.first_dyn = None 
        self.last_dyn = None 

        self.n_images = 0
        self.k_spring = 1

        self.elastic_band = [] 
        self.minim_options = None 

        self.minimizations = []
    
    def stetup_calculation(self, first_dyn, last_dyn, n_images, minim):
        """
        SETUP NEB
        =========

        This function setups a new neb for the calculation.
        It requires initial and final position of the elastic band, 
        the number of images as well as the number of images to be created.
        A linear interpolation will be done between the points.

        Note: distance is defined with atoms ordered. Therefore, be carefull, as
        changing the order of the atoms in the dynamical matrix can drastically change the result.

        Paramters
        ---------
            first_dyn : Phonons()
                The starting dynamical matrix
            last_dyn : Phonons()
                The end of the band.
            n_images : int
                The size of the elastic band (minimum 2, i.e. only starting and ending points)
            minim : SSCHA_Minimizer()
                The options for the sscha minimization
        """

        self.n_images = n_images
        self.first_dyn = first_dyn.Copy()
        self.last_dyn = last_dyn.Copy()
        self.elastic_band = []

        self.minim_options = minim

        # Interpolate in between
        for i in range(self.n_images):
            # Skip the first and the last
            if i == 0:
                self.elastic_band.append(self.first_dyn)
                continue 
            elif i == self.n_images - 1:
                self.elastic_band.append(self.last_dyn)
                continue 
            
            new_dyn = self.first_dyn.Copy()
            new_dyn.structure.coords += i * (self.last_dyn.structure.coords - self.first_dyn.structure.coords) / np.float64(n_images-1)
            for iq, q in enumerate(new_dyn.q_tot):
                new_dyn.dynmats[iq] += i * (self.last_dyn.dynmats[iq] - self.first_dyn.dynmats[iq]) / np.float64(n_images-1)

            self.elastic_band.append(new_dyn)


    def save_status(self, prefix):
        """
        SAVE THE NEB STATUS
        ===================

        This option will save the chain of dynamical matrices as

        PREFIX_neb_{:04d}_X

        where X is the Q point of the matrix, while the format string is the position
        in the neb chain.

        A directory PREFIX_ensemble will be created, that will store the last ensembles
        for all the chain, the population will be the index. The ensemble will be saved
        only if the minimizations have been initialized
        """

        for i in range(self.n_images):
            self.elastic_band[i].save_qe("{}_neb_{:04d}_".format(prefix, i))
        
        if len(self.minimizations) == self.n_images:
            os.makedirs("{}_ensemble".format(prefix))
            for i in range(self.n_images):
                self.minimizations[i].ensemble.save_bin("{}_ensemble".format(prefix), i)

        
    def run(self, n_pre_relax = 1):
        """
        RUN THE NEB
        ===========

        This function will run the NEB. 
        You can specify a thermalization option, to keep fixed the atomic
        position for some populations to optimize better the dynamical matrix.

        Paramters
        ---------
            n_pre_relax : int
                The  number of populations to relax only the dynamical matrices
                before actually starting the NEB
        """
        pass

    def NEB_CFG(self, grad_dyn, grad_struct):
        """
        The function that project the gradient of the structure into the 
        NEB one.
        """
        pass

    def single_iteration(self, use_neb = True):
        """
        Perform a single iteration over the whole band.
        
        Parameters
        ----------
            use_neb : bool
                If true the NEB minimization is performed,
                otherwise, atomic positions are kept fixed.
        """
        pass




            
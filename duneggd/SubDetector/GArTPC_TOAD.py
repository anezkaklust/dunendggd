""" GArTPC_TOAD.py

Author: Anežka Klustová, Imperial College London, a.klustova20@imperial.ac.uk

Based on GArTPC.py by J. Lopez, U. Colorado

A basic builder for a gas TPC for TOAD consisting of:
- a cylindrical chamber with two endcaps filled with gas
- a cylindrical pressure vessel with two endcaps 
- cathode
- field cage
- OROC

Gas region (cylindrical volume + endcaps) is an active volume.

TO DO:

- Offset of the terminator from the OROC?
- Beam window
- Make metal clamps?
- Include finer structure of the OROC - wedges?

"""

import gegede.builder
from duneggd.LocalTools import localtools as ltools
from duneggd.LocalTools import materialdefinition as materials
from gegede import Quantity as Q
from math import *


class GArTPCBuilder(gegede.builder.Builder):
    """ Class to build a gaseous argon TPC geometry with pressure vessel
        and inner components."

    Attributes:
        halfDimension: Dimensions of the entire bounding volume holding the geometry.
        ChamberRadius: Inner radius of the vacuum vessel.
        ChamberLength: Full length of the straight part of the vessel.
        GasType: Name of gas material.
        GasDensity: Density of gas (only used for custom gas mixes).
        Composition: Composition of a custom gas mixture.
        Drift: Drift axis (always z for a cylindrical TPC)
        SmallGap: A small distance to help prevent overlaps
    """

    def configure(self,chamberDimension,
                  halfDimension,GasType,drift='z',**kwargs):

        """ Set the configuration for the geometry.

            The keywords MaterialName and Density should only be used
            if Material is a dict-type rather than a string.

            Args:
                halfDimension: Half-dimensions of bounding volume.
                    Dict with keys 'rmin', 'rmax' and 'dz' (dz=half of length)
                chamberDimension: Outer dimensions of vauum vessel.
                    Dict. with keys 'r' and 'dz'
                GasType: Gas material. String if using a standard
                    material, dict in the form {material:mass_fraction,...}
                drift: The drift direction. (x, y, or z)
                kwargs: Additional keyword arguments. Allowed are:
                    TO DO - write list
        """

        self.halfDimension = halfDimension # bounding volume
        self.Material = 'Air' # material  of the boudning volume
        self.Drift = drift
        self.SmallGap = Q('0.001mm') # safety for overlaps

        # The gas
        if type(GasType)==str:
            # Set from a pre-defined material
            self.GasType = GasType
            self.GasDensity = None
            self.Composition = None
        else:
            # Set from a dictionary of materials & mass fractions
            comp = []
            for k in GasType:
                comp.append( (k,GasType[k]) )
            self.Composition = tuple(comp)
            self.GasType = kwargs['GasName']
            self.GasDensity = kwargs['Density']

        # Gaseous Chamber Details
        self.ChamberRadius = chamberDimension['r'] # inner radius
        self.ChamberLength = chamberDimension['dz'] # length of the straight part of the vessel
        self.RailLength = Q("625mm") # centered in the pressure vessel, shorted than the length of the vessel
        self.EndCapBulge = Q("32.5cm")

        # Pressure vessel
        self.pvThickness = Q("1cm")
        if 'pvThickness' in list(kwargs.keys()):
            self.pvThickness = kwargs['pvThickness']
        self.pvMaterial = "Steel"
        if 'pvMaterial' in list(kwargs.keys()):
            self.pvMaterial= kwargs['pvMaterial']

        # Offsets in z
        # Offset from the end of the rail on the OROC side
        self.TerminatorOffset = Q("284mm")
        if "TerminatorOffset" in list(kwargs.keys()):
            self.TerminatorOffset = kwargs['TerminatorOffset']
        # Distance between the terminator and the OROC
        self.TerminatorOROCoffset = Q("10mm")
        if "TerminatorOROCoffset" in list(kwargs.keys()):
           self.TerminatorOROCoffset = kwargs['TerminatorOROCoffset']
        # Distance between the terminator and the first ring of the field cage
        self.fcTerminatorFirstRingOffset = Q("43mm")
        if "fcTerminatorFirstRingOffset" in list(kwargs.keys()):
            self.cTerminatorFirstRingOffset = kwargs['fcTerminatorFirstRingOffset']

        # Field Cage (Rings)
        self.fcInnerRadius = Q("55.5cm")
        if "fcInnerRadius" in list(kwargs.keys()):
            self.fcInnerRadius= kwargs['fcInnerRadius']
        self.fcOuterRadius = Q("56.1cm")
        if "fcOuterRadius" in list(kwargs.keys()):
            self.fcOuterRadius = kwargs['fcOuterRadius']
        self.fcNRings = "12"
        if "fcNRings" in list(kwargs.keys()):
            self.fcNRings = kwargs['fcNRings']
        self.fcRingSpacing = Q("2.5cm")
        if "fcRingSpacing" in list(kwargs.keys()):
            self.fcRingSpacing = kwargs['fcRingSpacing']
        self.fcMaterial = "Steel"
        if "fcMaterial" in list(kwargs.keys()):
            self.fcMaterial = kwargs['fcMaterial']
        self.fcRingThickness =  Q("1cm")
        if "fcRingThickness" in list(kwargs.keys()):
            self.fcRingThickness = kwargs['fcRingThickness']

        # OROC Main trapezoid
        # x, y, z correpond to x, y, z orientation of the bounding volume
        # Upper, Lower, Thickness, Height correspond to the actual orientation after placing
        # it in the chamber (after rotation about 270 deg around x)
        self.orocMaterial=  "Aluminum"
        if "orocMaterial" in list(kwargs.keys()):
            self.orocMaterial = kwargs['orocMaterial']
        c = Q("467.747mm")
        if "orocUpper_dx1" in list(kwargs.keys()):
            self.orocUpper_dx1 = kwargs['orocUpper_dx1']
        self.orocLower_dx2 =  Q("870.478mm")
        if "orocLower_dx2" in list(kwargs.keys()):
            self.orocLower_dx2 = kwargs['orocLower_dx2']
        self.orocThickness_dy1 =  Q("21.8mm") 
        if "orocThickness_dy1" in list(kwargs.keys()):
            self.orocThickness_dy1 = kwargs['orocThickness_dy1']
        self.orocThickness_dy2 =  self.orocThickness_dy1 # same thickness everywhere
        if "orocThickness_dy2" in list(kwargs.keys()):
           self.orocThickness_dy2  = kwargs['orocThickness_dy2']
        self.orocHeight_dz =  Q("1142mm")
        if "orocHeight_dz" in list(kwargs.keys()):
            self.orocHeight_dz = kwargs['orocHeight_dz']
        self.orocOffset_y = Q("56mm") # calculated
        # OROC Main offset in y (calculated)
        if "orocOffset_y" in list(kwargs.keys()):
            self.orocOffset_y = kwargs['orocOffset_y']

        # OROC pad frame
        # same material, also Aluminum
        self.orocFrameThickness_dy1  =  Q("62mm")
        if "orocFrameThickness_dy1" in list(kwargs.keys()):
           self.orocFrameThickness_dy1 = kwargs['orocFrameThickness_dy1']
        self.orocFrameThickness_dy2  =  self.orocFrameThickness_dy1 
        if "orocFrameThickness_dy2" in list(kwargs.keys()):
           self.orocFrameThickness_dy2 = kwargs['orocFrameThickness_dy2']
        self.orocFrameHeight_dz = Q("1027mm")
        self.orocFrameVerticalDivisionThickness_dx1 = Q("16mm")
        if "orocFrameVerticalDivisionThickness_dx1" in list(kwargs.keys()):
           self.orocFrameVerticalDivisionThickness_dx1 = kwargs['orocFrameVerticalDivisionThickness_dx1']
        self.orocFrameVerticalDivisionThickness_dx2 = self.orocFrameVerticalDivisionThickness_dx1  # same
        if "orocFrameVerticalDivisionThickness_dx2" in list(kwargs.keys()):
           self.orocFrameVerticalDivisionThickness_dx2 = kwargs['orocFrameVerticalDivisionThickness_dx2']
        self.orocFrameVerticalDivisionDepth_dy1 = Q("42.5mm")
        if "orocFrameVerticalDivisionDepth_dy1" in list(kwargs.keys()):
           self.orocFrameVerticalDivisionDepth_dy1 = kwargs['orocFrameVerticalDivisionDepth_dy1']
        self.orocFrameVerticalDivisionDepth_dy2 = self.orocFrameVerticalDivisionDepth_dy1  # same
        if "orocFrameVerticalDivisionDepth_dy2" in list(kwargs.keys()):
           self.orocFrameVerticalDivisionDepth_dy2 = kwargs['orocFrameVerticalDivisionDepth_dy2']
        self.orocFrameOffset_y = Q("77.5mm")  # calculated
        # OROC Frame offset in y (calculated)
        if "orocFrameOffset_y " in list(kwargs.keys()):
            self.orocFrameOffset_y  = kwargs['orocFrameOffset_y ']

        # Cathode
        # Wire mesh ~ 25 microns thick
        self.CathodeMaterial = "Steel"
        if "CathodeMaterial" in list(kwargs.keys()):
            self.CathodeMaterial= kwargs['CathodeMaterial']
        self.CathodeThickness = Q("0.025mm") # thickness of the mesh is about 25 microns
        if "CathodeThickness" in list(kwargs.keys()):
            self.c_thickness= kwargs['CathodeThickness']
        self.CathodeInnerRadius =  Q("0mm")
        if "CathodeInnerRadius" in list(kwargs.keys()):
            self.c_rInner = kwargs['CathodeInnerRadius']
        self.CathodeOuterRadius =  Q("56cm") # solid ring
        if "CathodeOuterRadius" in list(kwargs.keys()):
            self.CathodeOuterRadius = kwargs['CathodeOuterRadius']

        # Cathode Holder
        self.CathodeHolderMaterial = "Steel"
        if "CathodeHolderMaterial" in list(kwargs.keys()):
            self.CathodeHolderMaterial = kwargs['CathodeHolderMaterial']
        self.CathodeHolderThickness = Q("3mm") #steel ring is 3mm wide
        if "CathodeHolderThickness" in list(kwargs.keys()):
            self.CathodeHolderThickness  = kwargs['CathodeHolderThickness']
        self.CathodeHolderInnerRadius = Q("56cm") # inner diameter is 112cm
        if "CathodeHolderInnerRadius" in list(kwargs.keys()):
            self.CathodeHolderInnerRadius = kwargs['CathodeHolderInnerRadius']
        self.CathodeHolderOuterRadius= Q("59cm") # outer diameter 118cm
        if "CathodeHolderOuterRadius" in list(kwargs.keys()):
            self.CathodeHolderOuterRadius = kwargs['CathodeHolderOuterRadius']

        # Terminator 
        # encloses the drif region on the OROC side
        self.TerminatorMaterial = "Steel"
        if "TerminatorMaterial" in list(kwargs.keys()):
            self.TerminatorMaterial = kwargs['TerminatorMaterial']
        self.TerminatorThickness = Q("1mm")
        if "TerminatorThickness" in list(kwargs.keys()):
            self.TerminatorThickness  = kwargs['TerminatorThickness']
        self.TerminatorInnerRadius = Q("0cm")
        if "TerminatorInnerRadius" in list(kwargs.keys()):
            self.TerminatorInnerRadius = kwargs['TerminatorInnerRadius']
        self.TerminatorOuterRadius = Q("56cm")
        if "TerminatorOuterRadius" in list(kwargs.keys()):
            self.TerminatorOuterRadius = kwargs['TerminatorOuterRadius']

        # Terminator Holder
        self.TerminatorHolderMaterial = "PVC"
        if "TerminatorHolderMaterial" in list(kwargs.keys()):
            self.TerminatorHolderMaterial = kwargs['TerminatorHolderMaterial']
        self.TerminatorHolderThickness = Q("1mm")
        if "TerminatorHolderThickness" in list(kwargs.keys()):
            self.TerminatorHolderThickness  = kwargs['TerminatorHolderThickness']
        self.TerminatorHolderInnerRadius = Q("56cm")
        if "TerminatorHolderInnerRadius" in list(kwargs.keys()):
            self.TerminatorHolderInnerRadius = kwargs['TerminatorHolderInnerRadius']
        self.TerminatorHolderOuterRadius = Q("59cm")
        if "TerminatorHolderOuterRadius" in list(kwargs.keys()):
            self.TerminatorHolderOuterRadius = kwargs['TerminatorHolderOuterRadius']



    def construct(self,geom):
        """ Construct the geometry.

        The standard geometry consists of a cylindrical vessel
        filled with gas with two endcaps enclosed in the pressure vessel.
        TOAD consists of cathode, field cage, OROC and the termiantor.

        args:
            geom: The geometry

        """

        # If using a custom gas, define here
        if self.Composition is not None:
            geom.matter.Mixture(self.GasType,
                                density=self.GasDensity,
                                components=self.Composition)

        main_lv, main_hDim = ltools.main_lv(self,geom, 'Tubs')
        print('GasTPCBuilder::construct()')
        print('main_lv = '+main_lv.name)
        self.add_volume(main_lv)
        print(ltools.getShapeDimensions(main_lv,geom))

        ###########################################
        ##### GASEOUS CHAMBER (STRAIGHT PART) #####
        ###########################################
        # Construct the chamber
        print("-------------------------------------------------")
        print("Construct gaseous chamber - barrel")
        tpc_chamber_shape = geom.shapes.Tubs('TPCChamber',
                                       rmax = self.ChamberRadius,
                                       dz = self.ChamberLength*0.5)
        tpc_chamber_lv = geom.structure.Volume('volTPCChamber',
                                               material=self.GasType,
                                               shape=tpc_chamber_shape)
        
        # The gas volumes are sensitive detectors
        tpc_chamber_lv.params.append(('SensDet',"Gas Barrel vol"))

        # Place into main LV
        pos = [Q('0m'),Q('0m'),Q('0m')]
        tpc_chamber_pos = geom.structure.Position('TPCChamber_pos',
                                                  pos[0],pos[1],pos[2])
        tpc_chamber_pla = geom.structure.Placement('TPCChamber_pla',
                                                   volume=tpc_chamber_lv,
                                                   pos=tpc_chamber_pos)#,
                                                   #rot = main_rot)
        main_lv.placements.append(tpc_chamber_pla.name)

        
        ###########################################
        ######## Pressure vessel Barrel ###########
        ###########################################
        print("-------------------------------------------------")
        print("Construct pressure vessel - barrel")

        pv_rInner = self.ChamberRadius
        pvHalfLength = self.ChamberLength/2
        pv_rmin = pv_rInner
        pv_rmax = pv_rmin + self.pvThickness

        # build the pressure vessel barrel
        pvb_name = "PV_Barrel"
        pvb_shape = geom.shapes.Tubs(pvb_name, rmin=pv_rmin, rmax=pv_rmax, dz=pvHalfLength, sphi="0deg", dphi="360deg")
        pvb_vol = geom.structure.Volume("vol"+pvb_name, shape=pvb_shape, material=self.pvMaterial)

        pvb_pla = geom.structure.Placement("PVBarrel"+"_pla", volume=pvb_vol)#, rot = main_rot)
        # Place it in the main lv
        main_lv.placements.append(pvb_pla.name)

        ###########################################
        ############ GASEOUS ENDCAPS  #############
        ###########################################
        print("-------------------------------------------------")
        print("Construct gaseous endcaps")
        # Position
        pv_rInner = self.ChamberRadius
        pv_rmin = sqrt((pv_rInner/Q("1mm"))**2)*Q("1mm")
        h = self.EndCapBulge
        x = pv_rmin
        q = ((h/Q("1mm"))**2 + (x/Q("1mm"))**2)
        R = q/(2*h/Q("1mm"))*Q("1mm")
        xpos = self.ChamberLength/2-(R-h)
        print("PV Endcap gas put at xpos", xpos)

        pv_rInner = self.ChamberRadius + self.SmallGap
        pvHalfLength = self.ChamberLength/2
        pv_rmin = pv_rInner

        h = self.EndCapBulge
        x = pv_rmin
        q = ((h/Q("1mm"))**2 + (x/Q("1mm"))**2)
        R = q/(2*h/Q("1mm"))*Q("1mm")
        dtheta = asin( 2*(h/Q("1mm"))*(x/Q("1mm"))/q)

        print("h, x, q, R, dtheta = ", h, x, q, R, dtheta)

        pvec_name =  "Endcap_Gas"
        pvec_shape = geom.shapes.Sphere(pvec_name + "_shape", rmin=Q("0cm"), rmax=R-self.SmallGap, sphi="0deg", dphi="360deg", stheta="0deg", dtheta=dtheta)

        # build the pressure vessel barrel
        pvb_name = "Gas_to_subtract"
        pvb_shape = geom.shapes.Tubs(pvb_name,rmax=self.ChamberRadius-self.pvThickness, dz=self.ChamberLength/2 + abs(xpos) + self.SmallGap)
        
        endcaps_shape = geom.shapes.Boolean('Gas_Endcaps_subtract_inner_barrel', type='subtraction', first=pvec_shape, second=pvb_shape)
        
        pvec_vol = geom.structure.Volume("Gas_Endcaps_vol", shape=endcaps_shape, material=self.GasType)
        
        # The gas volumes are sensitive detectors
        pvec_vol.params.append(('SensDet',"Gas_Endcaps_vol"))
        
        for side in ["L", "R"]:
            yrot = "0deg" if side == 'L' else "180deg"
            if side == 'R':
                xpos = -xpos
            # print("xpos = ", xpos)
            pvec_rot = geom.structure.Rotation("Gas_Endcap_"+side+"_rot", y=yrot)
            pvec_pos = geom.structure.Position("Gas_Endcap_"+side+"_pos", z=xpos)
            pvec_pla = geom.structure.Placement("Gas_Endcap_"+side+"_pla", volume=pvec_vol, pos=pvec_pos, rot=pvec_rot)
            main_lv.placements.append(pvec_pla.name)

        ###########################################
        ##### ENDCAPS VOLUME Pressure vessel ######
        ###########################################
        print("-------------------------------------------------")
        print("Construct pressure vessel endcaps")
        # Position
        pv_rInner = self.ChamberRadius
        pv_rmin = sqrt((pv_rInner/Q("1mm"))**2)*Q("1mm")
        h = self.EndCapBulge
        x = pv_rmin
        q = ((h/Q("1mm"))**2 + (x/Q("1mm"))**2)
        R = q/(2*h/Q("1mm"))*Q("1mm")
        xpos = self.ChamberLength/2-(R-h)
        print("PV Endcap put at xpos", xpos)

        # build the pressure vessel endcaps
        pv_rInner = self.ChamberRadius
        pvHalfLength = self.ChamberLength/2
        pv_rmin = pv_rInner + self.SmallGap
        pv_rmax = pv_rmin + self.pvThickness
        
        h = self.EndCapBulge
        x = pv_rmin
        q = ((h/Q("1mm"))**2 + (x/Q("1mm"))**2)
        R = q/(2*h/Q("1mm"))*Q("1mm")
        dtheta = asin( 2*(h/Q("1mm"))*(x/Q("1mm"))/q)
        
        print("h, x, q, R, dtheta = ", h, x, q, R, dtheta)

        pvec_name = "PV_Endcap"
        pvec_shape = geom.shapes.Sphere(pvec_name, rmin=R, rmax=R + self.pvThickness, sphi="0deg", dphi="360deg", stheta="0deg", dtheta=dtheta)
        pvec_vol = geom.structure.Volume("vol_"+pvec_name, shape=pvec_shape, material=self.pvMaterial)
        
        for side in ["L", "R"]:
            yrot = "0deg" if side == 'L' else "180deg"
            if side == 'R':
                xpos = -xpos
            # print("xpos = ", xpos)
            pvec_rot = geom.structure.Rotation("PV_Endcap_"+side+"_rot", y=yrot)
            pvec_pos = geom.structure.Position("PV_Endcap_"+side+"_pos", z=xpos)
            pvec_pla = geom.structure.Placement("PV_Endcap_"+side+"_pla", volume=pvec_vol, pos=pvec_pos, rot=pvec_rot)
            main_lv.placements.append(pvec_pla.name)

        print("-------------------------------------------------")
        print("Construct stuff inside the chamber/vessel")
        print("-------------------------------------------------")

        self.construct_tpcs(geom, tpc_chamber_lv)

    def construct_tpcs(self,geom,lv):
        """ Construct the two TPCs along with their
        field cages and readout plane

        args:
            geom: The geometry:
            tpc_chamber_lv: The vessel gas volume where
                want to stuff
        """

        pos1 = []
        pos2 = []
        rot1 = []
        rot2 = []
        rot3 = []
        rot4 = []

        if self.Drift == 'y':

            pos1 = [0,1,0]
            pos2 = [0,-1,0]
            rot1 = [Q('270deg'),Q('0deg'),Q('0deg')]
            rot2 = [Q('90deg'),Q('0deg'),Q('0deg')]

        elif self.Drift == 'x':

            pos1 = [1,0,0]
            pos2 = [-1,0,0]
            rot1 = [Q('0deg'),Q('90deg'),Q('0deg')]
            rot2 = [Q('270deg'),Q('deg'),Q('0deg')]
            rot3 = [Q('30deg'),Q('0deg'),Q('0deg')]

        else:
            pos0 = [0,0,0] # origin
            # x is left and right within the vessel
            # y is up and bottom within the vessel
            # z is towards and away within the vessel

            rot1 = [Q('0deg'),Q('0deg'),Q('0deg')]
            rot2 = [Q('270deg'),Q('0deg'),Q('0deg')]
            rot3 = [Q('0deg'),Q('0deg'),Q('-10deg')]
            rot4 = [Q('0deg'),Q('0deg'),Q('10deg')]


        # have to define that the distance in between the places
        # first ring offset 
        ringOffset = self.RailLength/2 - self.TerminatorOffset - self.TerminatorThickness - self.fcTerminatorFirstRingOffset - self.fcRingThickness/2
        print("First ring offset from the center (from the OROC side) = " + str(ringOffset))
        # Based on 43 mm offset between the terminator and the first ring

        self.construct_cathode(geom,"Cathode",pos0,rot1,lv)
        self.construct_fieldcagering(geom,"FC_ringN_" + str(0),pos0, ringOffset, rot1,lv)
        print("Ring " + str(0) + " position in z = " + str(ringOffset))
        
        for ring in range(1,int(self.fcNRings)):
            self.construct_fieldcagering(geom,"FC_ringN_" + str(ring),pos0, ringOffset - ring*self.fcRingSpacing, rot1,lv)
            print("Ring " + str(ring) + " position in z = " + str(ringOffset - ring*self.fcRingSpacing))
    
        self.construct_terminator_holder(geom,"Terminator_holder",pos0,rot1,lv)
        self.construct_terminator_bottom(geom,"Terminator_bottom",pos0,rot1,lv)
        self.construct_terminator_left(geom,"Terminator_left",pos0,rot3,lv) 
        self.construct_terminator_right(geom,"Terminator_right",pos0,rot4,lv)
        self.construct_oroc_main(geom,"OROC_alu_main",pos0,rot2,lv) 
        self.construct_oroc_frame(geom,"OROC_alu_frame",pos0,rot2,lv)

    def construct_cathode(self, geom, name, pos_vec,rot, lv):
        # Cathode is 10.5 mm from the end of the rails (on the cathode side)
        # Mesh thicknes is about 25 microns
        # Steel ring is 1 mm thick

        # Note we do not model the black delring that is 12.5 mm thick
        # which is immediately after the cathode (on the OROC side)
        # Details:
        # From the front end of the rails it is 601 mm to the cathode holder (black delrin)
        # Then the cathode holder is 12.5 mm, then the cathode
        # Then it is another 11. 5 mm to the far end of the rails (edited) 
        # Actually the cathode is attached to a steel ring 1mm width, so it is actually 614.5 mm from one end of the rails
        # And 10.5 mm from the other)

        '''Construct cathode.'''
        print("Construct cathode")

        offset_z = self.RailLength/2 - Q("10.5mm") - self.CathodeThickness/2
        print("Cathode position in z " + str(-offset_z))
        
        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1],pos[2] - offset_z)

        c_holder_rInner = self.CathodeHolderInnerRadius
        c_holder_rOuter = self.CathodeHolderOuterRadius
        c_holder_thickness = self.CathodeHolderThickness/2
        c_holder_material = self.CathodeHolderMaterial

        c_rInner = self.CathodeInnerRadius
        c_rOuter = self.CathodeOuterRadius
        c_material = self.CathodeMaterial
        c_dZthickness = self.CathodeThickness/2 #half width as usual

        print("1) Cathode Mesh Holder ")
        c_holder_shape = geom.shapes.Tubs(name+"cathode_holder_shape", rmax = c_holder_rOuter,rmin=c_holder_rInner, dz = c_holder_thickness)
        c_holder_vol = geom.structure.Volume(name+'cathode_holder_vol', shape=c_holder_shape, material=c_holder_material)
        c_holder_pla = geom.structure.Placement(name+'cathode_holder_place', volume=c_holder_vol, pos=tpc_pos, rot=tpc_rot)
        
        lv.placements.append(c_holder_pla.name)

        print("2) Cathode Mesh ")
        c_shape = geom.shapes.Tubs(name+"cathode_shape", rmax = c_rOuter,rmin=c_rInner, dz = c_dZthickness)
        c_vol = geom.structure.Volume(name+'cathode_vol', shape=c_shape, material=c_material)
        c_pla = geom.structure.Placement(name+'cathode_place', volume=c_vol, pos=tpc_pos, rot=tpc_rot)        
        
        lv.placements.append(c_pla.name)

    def construct_fieldcagering(self,geom,name,pos_vec, offset, rot,lv):
        '''Construct 1 field cage ring.'''
        print("Construct_Field_Cage " + name)
        
        # field cage position
        # field cage rotation
        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1],pos[2] + offset)
        
        fc_rInner = self.fcInnerRadius
        fc_rOuter = self.fcOuterRadius
        fc_dZring = self.fcRingThickness
        fc_material = self.fcMaterial # copper

        fc_shape = geom.shapes.Tubs(name+"field_cage_shape", rmax = fc_rOuter,rmin=fc_rInner, dz = fc_dZring)
        fc_vol = geom.structure.Volume(name+'field_cage_vol', shape=fc_shape, material=fc_material)
        fc_pla = geom.structure.Placement(name+'field_cage_placement', volume=fc_vol, pos=tpc_pos, rot=tpc_rot)
        
        lv.placements.append(fc_pla.name)
    
    def construct_terminator_holder(self, geom, name, pos_vec,rot, lv):
        '''Construct terminator holder'''
        print("Construct terminator holder")

        offset_z = self.RailLength/2 - self.TerminatorOffset - self.TerminatorThickness/2
        print("Terminator position in z " + str(offset_z))

        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1],pos[2]+offset_z)
        # off set readout plane by the thickness of the oroc trapezoid in negative z direction (inside of the vessel)
    
        t_holder_rInner = self.TerminatorHolderInnerRadius
        t_holder_rOuter = self.TerminatorHolderOuterRadius
        t_holder_thickness = self.TerminatorHolderThickness/2
        t_holder_material =  self.TerminatorHolderMaterial 
        # made of PVC

        t_holder_shape = geom.shapes.Tubs(name+"terminator_holder_shape", rmax = t_holder_rOuter,rmin=t_holder_rInner, dz = t_holder_thickness)
        t_holder_vol = geom.structure.Volume(name+'terminator_holder_vol', shape=t_holder_shape, material=t_holder_material)
        t_holder_pla = geom.structure.Placement(name+'terminator_holder_place', volume=t_holder_vol, pos=tpc_pos, rot=tpc_rot)
        
        lv.placements.append(t_holder_pla.name)

    def construct_terminator_bottom(self, geom, name, pos_vec,rot, lv):
        '''Construct bottom part of the terminator'''
        print("Construct terminator bottom")

        offset_z = self.RailLength/2 - self.TerminatorOffset - self.TerminatorThickness/2

        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1],pos[2]+offset_z)
        # off set readout plane by the thickness of the oroc trapezoid in negative z direction (inside of the vessel)
    
        t_rInner = self.TerminatorInnerRadius
        t_rOuter = self.TerminatorOuterRadius
        t_material = self.TerminatorMaterial 
        t_dZthickness = self.TerminatorThickness/2 #half width as usual

        # block 1
        oroc_dx1 = self.orocUpper_dx1 *2 # min (top)
        oroc_dx2 = self.orocUpper_dx1 *2  # still min (top)
        oroc_dy1 =  self.orocHeight_dz /2.2  # height
        oroc_dy2 =  self.orocHeight_dz /2.2  # height
        oroc_dz = self.TerminatorThickness/2 + Q("1cm")

        t_shape = geom.shapes.Tubs(name+"_terminator_shape_1", rmax = t_rOuter,rmin=t_rInner, dz = t_dZthickness, dphi = Q("180deg"), sphi = Q("180deg")) # lower circle
        oroc_shape = geom.shapes.Trapezoid(name + "_oroc_shape", dx1=oroc_dx1, dx2=oroc_dx2, dy1=oroc_dy1, dy2=oroc_dy2, dz=oroc_dz)

        # boolean
        subtract_terminator = geom.shapes.Boolean(name + '_terminator_subtract', type='subtraction', first=t_shape, second=oroc_shape)

        terminator_vol = geom.structure.Volume(name + "_terminator_vol_1", shape=subtract_terminator, material=t_material)
        terminator_pla = geom.structure.Placement(name+'_terminator_place_1', volume=terminator_vol, pos=tpc_pos, rot=tpc_rot)

        lv.placements.append(terminator_pla.name)

    def construct_terminator_left(self, geom, name, pos_vec,rot, lv):
        '''Construct left part of the terminator'''
        print("Construct terminator left")

        offset_z = self.RailLength/2 - self.TerminatorOffset - self.TerminatorThickness/2

        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1],pos[2]+offset_z)
        # off set readout plane by the thickness of the oroc trapezoid in negative z direction (inside of the vessel)
    
        t_rInner = self.TerminatorInnerRadius
        t_rOuter = self.TerminatorOuterRadius
        t_material =  self.TerminatorMaterial 
        t_dZthickness = self.TerminatorThickness/2 #half width as usual

        # block 1
        oroc_dx1 = self.orocUpper_dx1/1.35 # min (top)
        oroc_dx2 =  self.orocUpper_dx1/1.35 # still min (top)
        # change x to ensure the proper area is cut out
        oroc_dy1 =  self.orocHeight_dz/2  # height
        oroc_dy2 =  self.orocHeight_dz/2  # height
        oroc_dz = self.TerminatorThickness/2 + Q("1cm")

        t_shape = geom.shapes.Tubs(name+"_terminator_shape_2", rmax = t_rOuter,rmin=t_rInner, dz = t_dZthickness, dphi = Q("180deg"), sphi = Q("270deg")) # right circle
        oroc_shape = geom.shapes.Trapezoid(name + "_oroc_shape", dx1=oroc_dx1, dx2=oroc_dx2, dy1=oroc_dy1, dy2=oroc_dy2, dz=oroc_dz)

        # boolean
        subtract_terminator = geom.shapes.Boolean(name + '_terminator_subtract', type='subtraction', first=t_shape, second=oroc_shape)

        terminator_vol = geom.structure.Volume(name + "_terminato_vol_2", shape=subtract_terminator, material=t_material)
        terminator_pla = geom.structure.Placement(name+'_terminator_place_2', volume=terminator_vol, pos=tpc_pos, rot=tpc_rot)

        lv.placements.append(terminator_pla.name)

    def construct_terminator_right(self, geom, name, pos_vec,rot, lv):
        '''Construct right part of the terminator'''
        print("Construct terminator right")

        offset_z = self.RailLength/2 - self.TerminatorOffset - self.TerminatorThickness/2
        
        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1],pos[2]+offset_z)
        # off set readout plane by the thickness of the oroc trapezoid in negative z direction (inside of the vessel)
    
        t_rInner = self.TerminatorInnerRadius
        t_rOuter = self.TerminatorOuterRadius
        t_material =  self.TerminatorMaterial 
        t_dZthickness = self.TerminatorThickness/2 #half width as usual

        # block 1
        oroc_dx1 = self.orocUpper_dx1/1.35 # min (top)
        oroc_dx2 =  self.orocUpper_dx1/1.35 # still min (top)
        # change x to ensure the proper area is cut out
        oroc_dy1 =  self.orocHeight_dz/2  # height
        oroc_dy2 =  self.orocHeight_dz/2  # height
        oroc_dz = self.TerminatorThickness/2 + Q("1cm")

        t_shape = geom.shapes.Tubs(name+"terminator_shape_3", rmax = t_rOuter,rmin=t_rInner, dz = t_dZthickness, dphi = Q("180deg"), sphi = Q("90deg")) # left circle 
        oroc_shape = geom.shapes.Trapezoid(name + "oroc_shape", dx1=oroc_dx1, dx2=oroc_dx2, dy1=oroc_dy1, dy2=oroc_dy2, dz=oroc_dz)

        # boolean
        subtract_terminator = geom.shapes.Boolean(name + '_terminator_subtract', type='subtraction', first=t_shape, second=oroc_shape)

        terminator_vol = geom.structure.Volume(name + "_terminator_vol_3", shape=subtract_terminator, material=t_material)
        terminator_pla = geom.structure.Placement(name+'_terminator_place_3', volume=terminator_vol, pos=tpc_pos, rot=tpc_rot)

        lv.placements.append(terminator_pla.name)


    def construct_oroc_main(self, geom, name, pos_vec,rot, lv):
        '''Construct oroc main trapezoid'''
        print("Construct OROC main alu body")
        # top of OROC and vessel 76 mm
        # bottom of oroc and vessel 188 mm
        offset_y = self.orocOffset_y # calculated
        offset_z = self.RailLength/2 - self.TerminatorOffset + self.TerminatorOROCoffset + self.orocThickness_dy1/2
        print("OROC main position in y " + str(offset_y))
        print("OROC main position in z " + str(offset_z))

        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        #pos = [ x*Q("3.5cm") for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1] + offset_y ,pos[2] + offset_z)

        oroc_dx1 = self.orocUpper_dx1/2 # min (top)
        oroc_dx2 = self.orocLower_dx2/2 # max (bottom)
        oroc_dy1 = self.orocThickness_dy1/2 # thickness
        oroc_dy2 = self.orocThickness_dy2/2 # thickness
        oroc_dz = self.orocHeight_dz /2 # height
        oroc_material = self.orocMaterial
        
        oroc_shape = geom.shapes.Trapezoid(name + "_OROC_main_shape", dx1=oroc_dx1, dx2=oroc_dx2, dy1=oroc_dy1, dy2=oroc_dy2, dz=oroc_dz)
        oroc_vol = geom.structure.Volume(name + "_OROC_main_vol", shape=oroc_shape, material=oroc_material)
        oroc_pla = geom.structure.Placement(name+'_OROC_main_place', volume=oroc_vol, pos=tpc_pos, rot=tpc_rot)

        lv.placements.append(oroc_pla.name)

    def construct_oroc_frame(self, geom, name, pos_vec,rot, lv):
        '''Construct oroc frame alu trapezoid'''
        print("Construct OROC frame alu body")
        # top of OROC and vessel 76 mm
        # bottom of oroc and vessel 188 mm

        offset_y = self.orocFrameOffset_y # calculated
        offset_z = self.RailLength/2 - self.TerminatorOffset + self.TerminatorOROCoffset + self.orocThickness_dy1 + self.orocFrameThickness_dy1/2
        print("OROC frame position in y " + str(offset_y))
        print("OROC frame position in z " + str(offset_z))

        tpc_rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        pos = [ x*self.SmallGap for x in pos_vec]
        #pos = [ x*Q("3.5cm") for x in pos_vec]
        tpc_pos = geom.structure.Position(name+'_pos',pos[0],pos[1] + offset_y ,pos[2] + offset_z)

        # Trapezoid which is by 16.5mm smaller than main
        # 16.5 from technical drawings
        oroc_dx1 = self.orocUpper_dx1/2 - Q("16.5mm") # min (top)
        oroc_dx2 = self.orocLower_dx2/2 - Q("16.5mm") # max (bottom)
        oroc_dy1 = self.orocFrameThickness_dy1/2 # thickness
        oroc_dy2 = self.orocFrameThickness_dy2/2# thickness
        oroc_dz = self.orocFrameHeight_dz/2 # height
        oroc_material = self.orocMaterial

        # Measurements to subtract the inner side and leave just the fame
        # 20 mm thicknes of the frame from technical drawings
        oroc_dx1_2 = self.orocUpper_dx1/2 - Q("16.5mm") - Q("20mm") # min (top)
        oroc_dx2_2 = self.orocLower_dx2/2 - Q("16.5mm") - Q("20mm") # max (bottom)
        oroc_dy1_2 = self.orocFrameThickness_dy1/2 + self.SmallGap # thickness
        oroc_dy2_2 = self.orocFrameThickness_dy2/2 + self.SmallGap # thickness
        oroc_dz_2 = self.orocFrameHeight_dz/2 - Q("20mm") # height

        # vertical division
        ver_div_dx1 =  self.orocFrameVerticalDivisionThickness_dx1/2 # min (top)
        ver_div_dx2 =  self.orocFrameVerticalDivisionThickness_dx2/2# max (bottom)
        # thickness of the vertical division 16 mm
        # dimension/2 
        ver_div_dy1 = self.orocFrameVerticalDivisionDepth_dy1/2# thickness
        ver_div_dy2 = self.orocFrameVerticalDivisionDepth_dy2/2 # thickness
        ver_div_dz = oroc_dz_2 # height

        # Central vertical division
        vertical_division = geom.shapes.Trapezoid(name + "_vertical_division", dx1=ver_div_dx1, dx2=ver_div_dx2, dy1=ver_div_dy1, dy2= ver_div_dy2, dz= ver_div_dz)
        
        oroc_shape_outer = geom.shapes.Trapezoid(name + "_OROC_frame_shape_outer", dx1=oroc_dx1, dx2=oroc_dx2, dy1=oroc_dy1, dy2=oroc_dy2, dz=oroc_dz)
        oroc_shape_inner = geom.shapes.Trapezoid(name + "_OROC_frame_shape_inner", dx1=oroc_dx1_2, dx2=oroc_dx2_2, dy1=oroc_dy1_2, dy2=oroc_dy2_2, dz=oroc_dz_2)
        oroc_shape_1 = geom.shapes.Boolean(name + '_OROC_frame_shape_1', type='subtraction', first=oroc_shape_outer, second=oroc_shape_inner)
        oroc_shape = geom.shapes.Boolean(name + '_OROC_frame_shape', type='union', first=oroc_shape_1, second=vertical_division)
        oroc_vol = geom.structure.Volume(name + "_OROC_frame_vol", shape=oroc_shape, material=oroc_material)
        oroc_pla = geom.structure.Placement(name+'_OROC_frame_place', volume=oroc_vol, pos=tpc_pos, rot=tpc_rot)

        lv.placements.append(oroc_pla.name)
    

    '''
    # NOT USED
    def construct_oroc_back(self,geom,name,pos_vec,rot,lv):
        # top of OROC and vessel 76 mm
        # bottom of oroc and vessel 188 mm

        """ Construct the back of the OROC with the horizontal and vertical divisions."

        """ 
        # outer frame/pad frame thickness
        thickness = Q("62mm")/2
        offset = Q("62mm")/2 - self.oroc_dy1_thickness/2 + Q("10mm")
        #self.oroc_dy1_thickness/2 # safe gap
        rot = geom.structure.Rotation(name+'_rot',rot[0],rot[1],rot[2])
        #pos = [ x*Q("3.5cm") for x in pos_vec]
        pos = [ x*self.SmallGap for x in pos_vec]
        # central
        pos0 = geom.structure.Position(name+'_pos0',pos[0],pos[1] + Q("56mm") ,pos[2] + offset)
        # horizontal division upper
        pos1 = geom.structure.Position(name+'_pos1',pos[0],pos[1] + Q("56mm") + Q("265.167mm") ,pos[2] + offset)
        # horizontal division lower
        pos2 = geom.structure.Position(name+'_pos2',pos[0],pos[1] + Q("56mm") - Q("265.167mm"),pos[2] + offset)
        # smaller horizontal divisions (x4)
        pos3 = geom.structure.Position(name+'_pos3',pos[0],pos[1] + Q("56mm") + Q("265.167mm") + Q("265.167mm")/2,pos[2] + offset)
        pos4 = geom.structure.Position(name+'_pos4',pos[0],pos[1] + Q("56mm") + Q("265.167mm")/2,pos[2] + offset)
        pos5 = geom.structure.Position(name+'_pos5',pos[0],pos[1] + Q("56mm") - Q("265.167mm")/2,pos[2] + offset)
        pos6 = geom.structure.Position(name+'_pos6',pos[0],pos[1] + Q("56mm") - Q("265.167mm") - Q("265.167mm")/2 ,pos[2] + offset)
        # off set readout plane by the thickness of the oroc trapezoid in positive z direction (outside of the vessel)

        # outer frame
        oroc_dx1 = self.oroc_dx1_upper/2 # min (top)
        oroc_dx2 = self.oroc_dx2_lower/2 # max (bottom)
        #oroc_dy1 = self.oroc_dy1_thickness/2 # thickness
        #oroc_dy2 = self.oroc_dy2_thickness/2 # thickness
        oroc_dy1 = thickness/2 # thickness
        oroc_dy2 = thickness/2 # thickness
        # how thick is the frame???
        oroc_dz = self.oroc_dz_height/2 # height
        oroc_material = self.oroc_material

        # inner frame
        oroc_dx1_2 = oroc_dx1  - Q("5cm") # min (top)
        oroc_dx2_2 = oroc_dx2 - Q("5cm")# max (bottom)
        oroc_dy1_2 = oroc_dy1  + Q("0.5cm") # thickness
        oroc_dy2_2 = oroc_dy1 + Q("0.5cm") # thickness
        # thickness must be bigger than outer so we get empty space inside for the subtraction
        oroc_dz_2 = Q("1008.99mm")/2
        oroc_dz - Q("5cm") # height
        # what is the height of the inner region
        
        # vertical division
        ver_div_dx1 =  Q("16mm")/2 # min (top)
        ver_div_dx2 =  Q("16mm")/2# max (bottom)
        # what is the thicknes of the vertical division? it is 16 mm
        # dimension/2 
        ver_div_dy1 = oroc_dy1 # thickness
        ver_div_dy2 = oroc_dy2 # thickness
        ver_div_dz = oroc_dz_2 # height

        # central horizontal division 1
        hor_div_dx1 =  Q("552.64mm")/2 # Q("582.64mm")/2
        # top dimension (min)
        hor_div_dx2 =  Q("560.1mm")/2 #Q("590.1mm")/2 # max (bottom)
        # dimension/2 
        hor_div_dy1 = self.oroc_dy1_thickness/2 # thickness
        hor_div_dy2 = self.oroc_dy2_thickness/2 # thickness


        hor_div_dy1 = thickness/2 # thickness
        hor_div_dy2 = thickness/2 # thickness
        hor_div_dz = Q("16mm")/2 # height

        # -----------
        # upper horizontal division 1
        hor_div_dx1_1 =  Q("451.601mm")/2#Q("491.601mm")/2
        # top dimension (min)
        hor_div_dx2_1 =  Q("457.453mm")/2 #Q("497.453mm")/2# max (bottom)

        # lower horizontal division 2
        hor_div_dx1_2 =  Q("668.298mm")/2#Q("678.298mm")/2
        # top dimension (min)
        hor_div_dx2_2 =  Q("674.162mm")/2#Q("684.162mm")/2# max (bottom)

        # -------------
        # smaller horizontal divisions (4 of them)
        hor_div_dz_2 = Q("4mm")/2 # height
        # most upper one
        hor_div_dx1_3 =  Q("398.185mm")/2 #Q("448.185mm")/2  # top dimension (min)
        hor_div_dx2_3 =  Q("399.724mm")/2#Q("449.724mm")/2 # max (bottom)

        # 2nd lower
        hor_div_dx1_4 =  Q("499.268mm")/2 #Q("539.268mm")/2  # top dimension (min)
        hor_div_dx2_4 =  Q("501.078mm")/2 # Q("541.078mm")/2 # max (bottom)

        # 3rd lower
        hor_div_dx1_5 =  Q("596.579mm")/2  #Q("636.579mm")/2  # top dimension (min)
        hor_div_dx2_5  =  Q("597.939mm")/2  #Q("637.939mm")/2 # max (bottom)

        # 4th lower (lowest)
        hor_div_dx1_6 =  Q("695.099mm")/2 #Q("725.099mm")/2  # top dimension (min)
        hor_div_dx2_6 =  Q("696.473mm")/2 # Q("726.473mm")/2 # max (bottom)


        shape_out = geom.shapes.Trapezoid(name + "oroc shape OUT ", dx1=oroc_dx1, dx2=oroc_dx2, dy1=oroc_dy1, dy2=oroc_dy2, dz=oroc_dz)
        shape_in = geom.shapes.Trapezoid(name + "oroc shape IN ", dx1=oroc_dx1_2, dx2=oroc_dx2_2, dy1=oroc_dy1_2, dy2= oroc_dy2_2, dz= oroc_dz_2)
        # boolean operations
        readout_plane = geom.shapes.Boolean(name + 'Readout plane shape ', type='subtraction', first=shape_out, second=shape_in)

        readout_plane_vol = geom.structure.Volume(name + "Readout plane vol ", shape=readout_plane, material=oroc_material)
        readout_plane_pla = geom.structure.Placement(name+"Readout plane place ", volume=readout_plane_vol, pos=pos0, rot=rot)
        lv.placements.append(readout_plane_pla.name)

        # central divisions
        # central horizontal and vertical divisions
        vertical_division = geom.shapes.Trapezoid(name + "vertical division ", dx1=ver_div_dx1, dx2=ver_div_dx2, dy1=ver_div_dy1, dy2= ver_div_dy2, dz= ver_div_dz)
        vertical_division_vol = geom.structure.Volume(name + "Vertical division volume ", shape=vertical_division, material=oroc_material)
        vertical_division_pla = geom.structure.Placement(name+"Vertical division place ", volume=vertical_division_vol, pos=pos0, rot=rot)
        lv.placements.append(vertical_division_pla.name)
        
        horizontal_division = geom.shapes.Trapezoid(name + "horizontal division ", dx1=hor_div_dx1, dx2=hor_div_dx2, dy1=hor_div_dy1, dy2= hor_div_dy2, dz= hor_div_dz)
        horizontal_division_vol = geom.structure.Volume(name + "Horizontal division volume ", shape=horizontal_division, material=oroc_material)
        horizontal_division_pla = geom.structure.Placement(name+"Horizontal division place ", volume=horizontal_division_vol, pos=pos0, rot=rot)
        lv.placements.append(horizontal_division_pla.name)

        # thicker horizontal divisions
        horizontal_division_upper = geom.shapes.Trapezoid(name + "horizontal division uppe r", dx1=hor_div_dx1_1, dx2=hor_div_dx2_1, dy1=hor_div_dy1, dy2= hor_div_dy2, dz= hor_div_dz)
        horizontal_division_upper_vol = geom.structure.Volume(name + "Horizontal division upper volume ", shape=horizontal_division_upper, material=oroc_material)
        horizontal_division_upper_pla = geom.structure.Placement(name+"Horizontal division upper place ", volume=horizontal_division_upper_vol, pos=pos1, rot=rot)
        lv.placements.append(horizontal_division_upper_pla.name)

        horizontal_division_lower = geom.shapes.Trapezoid(name + "horizontal division lower ", dx1=hor_div_dx1_2, dx2=hor_div_dx2_2, dy1=hor_div_dy1, dy2= hor_div_dy2, dz= hor_div_dz)
        horizontal_division_lower_vol = geom.structure.Volume(name + "Horizontal division lower volume ", shape=horizontal_division_lower, material=oroc_material)
        horizontal_division_lower_pla = geom.structure.Placement(name+"Horizontal division lower place ", volume=horizontal_division_lower_vol, pos=pos2, rot=rot)
        lv.placements.append(horizontal_division_lower_pla.name)

        # thinner horizontal divisions
        hor_div_1 = geom.shapes.Trapezoid(name + "Hor div 1 ", dx1=hor_div_dx1_3, dx2=hor_div_dx2_3, dy1=hor_div_dy1, dy2= hor_div_dy2, dz= hor_div_dz_2)
        hor_div_1_vol = geom.structure.Volume(name + "Hor div 1 volume ", shape=hor_div_1, material=oroc_material)
        hor_div_1_pla = geom.structure.Placement(name+"Hor div 1 place ", volume=hor_div_1_vol, pos=pos3, rot=rot)
        lv.placements.append(hor_div_1_pla.name)

        hor_div_2 = geom.shapes.Trapezoid(name + "Hor div 2 ", dx1=hor_div_dx1_4, dx2=hor_div_dx2_4, dy1=hor_div_dy1, dy2= hor_div_dy2, dz= hor_div_dz_2)
        hor_div_2_vol = geom.structure.Volume(name + "Hor div 2 volume ", shape=hor_div_2, material=oroc_material)
        hor_div_2_pla = geom.structure.Placement(name+"Hor div 2 place ", volume=hor_div_2_vol, pos=pos4, rot=rot)
        lv.placements.append(hor_div_2_pla.name)

        hor_div_3 = geom.shapes.Trapezoid(name + "Hor div 3 ", dx1=hor_div_dx1_5, dx2=hor_div_dx2_5, dy1=hor_div_dy1, dy2= hor_div_dy2, dz= hor_div_dz_2)
        hor_div_3_vol = geom.structure.Volume(name + "Hor div 3 volume ", shape=hor_div_3, material=oroc_material)
        hor_div_3_pla = geom.structure.Placement(name+"Hor div 3 place ", volume=hor_div_3_vol, pos=pos5, rot=rot)
        lv.placements.append(hor_div_3_pla.name)

        hor_div_4 = geom.shapes.Trapezoid(name + "Hor div 4 ", dx1=hor_div_dx1_6, dx2=hor_div_dx2_6, dy1=hor_div_dy1, dy2= hor_div_dy2, dz= hor_div_dz_2)
        hor_div_4_vol = geom.structure.Volume(name + "Hor div 4 volume ", shape=hor_div_4, material=oroc_material)
        hor_div_4_pla = geom.structure.Placement(name+"Hor div 4 place ", volume=hor_div_4_vol, pos=pos6, rot=rot)
        lv.placements.append(hor_div_4_pla.name)
    '''
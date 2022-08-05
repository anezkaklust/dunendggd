#!/usr/bin/env python
'''
TOAD_Builder: Builds the multi purpose tracker
'''

import gegede.builder
from duneggd.LocalTools import localtools as ltools
from duneggd.LocalTools import materialdefinition as materials
from gegede import Quantity as Q
from math import asin, sqrt

class TOAD_Builder(gegede.builder.Builder):
    '''
    Build a concept of the TOAD detector. This class directly
    sub-builders for the GArTPC and the Pressure Vessel.

    Arguments:
    buildGarTPC: Flag to build the GArTPC
    buildPV: Flag to build the Pressure Vessel

    '''

    defaults=dict( innerBField="0.0 T, 0.0 T, 0.0 T",
                   TPCStepLimit = "10 mm",
                   buildGarTPC=True,
                   space=Q("10cm")
                   )

    #^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^
    def construct(self, geom):

        ''' Top level volume (MPD) - It is rotated later in the cavern (x, y, z) -> (z, y, x)'''
        
        #  to define materials (!!)
        materials.define_materials(geom)

        dx_main=Q('8m')+self.space #dimension along the beam
        dy_main=Q('8m')+self.space #dimension in height
        dz_main=Q('8m')+self.space #dimension perp to the beam

        print("Dimension of the MPD in along the beam ", dx_main*2, " dimension in height ", dy_main*2, " and dimension perp to the beam ", dz_main*2)

        main_shape = geom.shapes.Box('MPD', dx=dx_main, dy=dy_main, dz=dz_main)
        main_lv = geom.structure.Volume('vol'+main_shape.name, material='Air', shape=main_shape)

        self.add_volume(main_lv)

        ##### build a fake volume that contains the MPD without the magnet to define the magnetized volume correctly ###
        # fake_lv = self.buildMagnetizedVolume(main_lv, geom)

        ######### build the TPC       ##########################
        # use GArTPC Builder, but "disable" the cyrostat by tweaking
        # EndcapThickness, WallThickness, and ChamberMaterial
        # do that in the cfg file
        if self.buildGarTPC:
            print("Adding TPC to main volume")
            self.build_gartpc(main_lv, geom)
        
        # x is the dimension along the beam
        # y is the height
        # z is the dimension along the electric field

        # garsoft input is
        # x along the electric field
        # z along the beam
        # rotate X to Z
        #self.RockRotation()

        ######### build the pressure vessel  ###################
        # Build the Pressure Vessel using the parameters in the cfg file
        # The PV consists of a cylinder for the Barrel and
        # the intersection of the cylinder and a sphere for the Endcaps
        #if self.buildPV:
        #    print("Adding PV to main volume")
        #    self.build_pressure_vessel(main_lv, geom)

    
        return

    
    def build_gartpc(self, main_lv, geom):

        #Build TPC
        tpc_builder = self.get_builder('GArTPC')
        if tpc_builder == None:
            return

        tpc_vol = tpc_builder.get_volume()
        print(tpc_vol)
        # Add the magnetic field to the volume
        #tpc_vol.params.append(("BField", self.innerBField))
        #tpc_vol.params.append(("StepLimit", self.TPCStepLimit))

        tpc_pla = geom.structure.Placement("GArTPC"+"_pla", volume=tpc_vol)
        # Place it in the main lv
        main_lv.placements.append(tpc_pla.name)

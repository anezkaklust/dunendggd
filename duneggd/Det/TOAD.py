#!/usr/bin/env python
'''
TOAD_Builder: Builds the The Teststand of an Overpressure Argon Detector (TOAD) geometry
'''

import gegede.builder
from duneggd.LocalTools import localtools as ltools
from duneggd.LocalTools import materialdefinition as materials
from gegede import Quantity as Q
from math import asin, sqrt

class TOAD_Builder(gegede.builder.Builder):
    '''
    Build a concept of the TOAD detector - only one subclass.

    Arguments:
    buildGarTPC: Flag to build the GArTPC
    '''

    defaults=dict( buildGarTPC=True,
                   space=Q("10cm")
                   )

    #^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^
    def construct(self, geom):

        ''' Top level volume - It is rotated later in the cavern (x, y, z) -> (z, y, x)'''
        
        #  to define materials (!!)
        materials.define_materials(geom)

        dx_main=Q('8m')+self.space #dimension along the beam
        dy_main=Q('8m')+self.space #dimension in height
        dz_main=Q('8m')+self.space #dimension perp to the beam

        print("Dimension of the top level volume in along the beam ", dx_main*2, " dimension in height ", dy_main*2, " and dimension perp to the beam ", dz_main*2)

        main_shape = geom.shapes.Box('MPD', dx=dx_main, dy=dy_main, dz=dz_main )
        main_lv = geom.structure.Volume('vol'+main_shape.name, material='Air', shape=main_shape)

        self.add_volume(main_lv)


        ######### build the TOAD TPC + pressure vessel ##########################
        # use GArTPC TOAD Builder
        # define dimensions and materials in the cfg file
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

        return

    
    def build_gartpc(self, main_lv, geom):

        #Build TPC
        tpc_builder = self.get_builder('GArTPC')
        if tpc_builder == None:
            return

        tpc_vol = tpc_builder.get_volume()
        print(tpc_vol)

        rot = geom.structure.Rotation("Rotation (x<->z)", y="90deg")

        tpc_pla = geom.structure.Placement("GArTPC"+"_pla", volume=tpc_vol, rot = rot )
        # Place it in the main lv
        main_lv.placements.append(tpc_pla.name)

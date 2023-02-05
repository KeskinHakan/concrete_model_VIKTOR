# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:37:14 2022

@author: hakan
"""

from viktor.parametrization import Parametrization, Section, NumberField, OptionField,Tab,TextField, DownloadButton, IsEqual, Lookup, And, IsFalse,BooleanField

class LogoParametrization(Parametrization):
    tab_1 = Tab('Analysis Inputs')
    tab_1.geometry = Section('Geometry')
    tab_1.geometry.sectiontype= OptionField("Section Analysis", options=["Existing", "Retrofit"], default="Retrofit")
    tab_1.geometry.width = NumberField('Width', suffix='mm', default=500)
    tab_1.geometry.depth = NumberField('Depth', suffix='mm', default=500)
    tab_1.geometry.cover = NumberField('Concrete Cover', suffix='mm', default=35)
    tab_1.geometry.length = NumberField('Column Length', suffix='mm', default=5000)
    tab_1.geometry.degree = OptionField('Degree', options=[0, 90], default=0, description="Design direction of the system - 0: X Direction / 90: Y Direction")
    
    tab_1.material = Section('Material')
    tab_1.material.steel = OptionField("Rebar Strength", suffix='MPa', options=[220, 420], default=420)
    tab_1.material.concrete = OptionField("Concrete Strength", suffix='MPa', options=[16, 20, 25, 30, 35, 40, 50], default=35)

    tab_1.material.coefficient = OptionField("Material Coefficient", options=["Nominal", "Expected"], default="Nominal")
    tab_1.material.gamma_concrete = NumberField("γconcrete",default=1)
    tab_1.material.gamma_steel = NumberField("γsteel",default=1)
    
    tab_1.longitudinal = Section('Longitudinal Rebar')     
    tab_1.longitudinal.diameter = NumberField('Diameter', suffix='mm', default=20, description="Longitudinal Rebar Diameter")
    tab_1.longitudinal.nx = NumberField('Number of Rebar Layer - X', suffix='-', default=3, description="Total rebar layer along to X Direction")
    tab_1.longitudinal.ny = NumberField('Number of Rebar Layer - Y', suffix='-', default=3, description="Total rebar layer along to Y Direction")
    
    tab_1.transverse = Section('Transverse Rebar')
    tab_1.transverse.confinementtype = OptionField("Confinement Type", options=["Unconfined", "Confined"], default="Confined")
    tab_1.transverse.diameter = NumberField('Diameter', suffix='mm', default=8, description="Transverse Rebar Diameter")
    tab_1.transverse.nlegx = NumberField('Number of Rebar - X', suffix='-', default=2, description="Total leg of the stirrups along to X Direction of the section")
    tab_1.transverse.nlegy = NumberField('Number of Rebar- Y', suffix='-', default=2, description="Total leg of the stirrups along to Y Direction of the section")
    tab_1.transverse.space = NumberField('s ', suffix='mm', default=120, description="Spacing of the Stirrups")
    tab_1.force = Section('Axial Force')
    tab_1.force.axialforce = NumberField('Axial Force', suffix='kN', default=500)
    
    tab_1.shearforce = Section('Shear Force')
    tab_1.shearforce.reducefactor = NumberField('D',default=2, description="Strength Reduced Factor")
    tab_1.shearforce.earthquake = NumberField('Shear Force (E)', suffix='kN', default=200, description="Earthquake Shear Force")
    tab_1.shearforce.combination = NumberField('Shear Force (G+Q+E)', suffix='kN', default=20, description="Gravity + Earthquake Shear Force")
    tab_1.shearforce.combination_d = NumberField('Shear Force (G+Q+D*E)', suffix='kN', default=50, description="Gravity + Increased Earthquake Shear Force with D")
    tab_1.shearforce.moment_bottom = NumberField('Moment Capacity (Bottom Joint)', suffix='kNm', default=350)
    tab_1.shearforce.moment_top = NumberField('Moment Capacity (Top Joint)', suffix='kNm', default=450)
    
    tab_1.shearforce.download_existing = DownloadButton('Download Moment Curvature (Existing)', method='download_mass_breakdown_existing')
    tab_1.shearforce.download_retrofit = DownloadButton('Download Moment Curvature (Retrofit)', method='download_mass_breakdown_retrofit')
    
    tab_2 = Tab('Concrete Analysis Inputs')
    tab_2.concreteproperties = Section('Distance Between Rebars')
    tab_2.concreteproperties.ax = NumberField('ax', suffix='mm', default=150)
    tab_2.concreteproperties.ay = NumberField('ay', suffix='mm', default=150)
    
    tab_3 = Tab('FRP Inputs')
    tab_3.frp = Section('Shear Resistance')
    # tab_3.frp._param_z_visible = And(
    # IsFalse(Lookup('tab_3.frp.param_x')),
    # IsEqual(Lookup('tab_3.frp.param_y'), 5)
    # )
    # tab_3.frp.param_x = BooleanField('X')
    # tab_3.frp.param_y = NumberField('Y')
    # tab_3.frp.param_z = NumberField('Z', visible=tab_3.frp._param_z_visible)
    
    tab_3.frp._sf_visible = IsEqual(Lookup('tab_3.frp.striptype'), "Uncontinuous")
    tab_3.frp._wf_visible = IsEqual(Lookup('tab_3.frp.striptype'), "Uncontinuous")
    
    tab_3.frp.striptype = OptionField("Strip Type", options=["Continuous", "Uncontinuous"], default="Continuous")

    tab_3.frp.numberoffrp = NumberField('nf', default=1, description="Number of Wrap")
    tab_3.frp.thicknessofwrap = NumberField('tf', suffix='mm',default=0.111, step = 0.001, num_decimals=3, description="Thickness") 
    tab_3.frp.widthofwrap = NumberField('wf', suffix='mm', default=500, visible=tab_3.frp._wf_visible, description="Width")
    tab_3.frp.youngmodulusofwrap = NumberField('Ef', suffix='N/mm2', default=230000, description="Young' Modulus")
    tab_3.frp.strainofwrap = NumberField('ef', suffix='-', default=0.004)
    tab_3.frp.ultimatestrainofwrap = NumberField('eu', suffix='-', default=0.021, description="Ultimate Strain")
    tab_3.frp.effectivedepth = NumberField('d', suffix='mm', default=500, description="Effective Depth")
    tab_3.frp.distanceofwraps = NumberField('sf', suffix='mm', default=500, visible=tab_3.frp._sf_visible, description="Distance Between FRPs")


    tab_3.axial_frp = Section('Axial and Flexural Resistance')
    tab_3.axial_frp.striptype = OptionField("Strip Type", options=["Continuous", "Uncontinuous"], default="Continuous")
    tab_3.axial_frp.sectiontype = OptionField("Section Type", options=["Circular Section", "Ellipse Section", "Rectangular Section"], default="Rectangular Section")
    tab_3.axial_frp.rc = NumberField('Rounding Radius - rc', suffix='mm', default=25)
    tab_3.axial_frp.numberoffrp = NumberField('nf', default=2, description="Number of Wrap")
    tab_3.axial_frp.thicknessofwrap = NumberField('tf', suffix='mm',default=0.33,step = 0.01, num_decimals=2, description="Thickness") 
    tab_3.axial_frp.widthofwrap = NumberField('wf', suffix='mm', default=150, visible=tab_3.frp._wf_visible, description="Width")
    tab_3.axial_frp.youngmodulusofwrap = NumberField('Ef', suffix='N/mm2', default=227500, description="Young' Modulus")
    tab_3.axial_frp.strainofwrap = NumberField('ef', suffix='-', default=0.004)
    tab_3.axial_frp.ultimatestrainofwrap = NumberField('eu', suffix='-', default=0.012, description="Ultimate Strain")
    tab_3.axial_frp.effectivedepth = NumberField('d', suffix='mm', default=500, description="Effective Depth")
    tab_3.axial_frp.distanceofwraps = NumberField('sf', suffix='mm', default=150, visible=tab_3.frp._sf_visible, description="Distance Between FRPs")
    
    # tab_3.axial_frp.pf = NumberField('Volumetric Ratio - pf', suffix='mm', default=0.05)

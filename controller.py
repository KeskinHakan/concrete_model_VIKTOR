from math import cos, sin, pi, sqrt, ceil
import pandas as pd
from viktor.core import ViktorController
from viktor.geometry import Line, Point, RectangularExtrusion, Material, Point, CartesianAxes, SquareBeam
from viktor.views import GeometryResult, GeometryView, DataResult, DataGroup, DataView,DataView, DataItem, PlotlyView, PlotlyResult, PlotlyAndDataView, PNGResult  
from .parametrization import LogoParametrization, DownloadButton
import matplotlib.pyplot as plt
from io import StringIO
from viktor.views import SVGView, SVGResult, PNGView
from openseespy.opensees import *
import openseespy.postprocessing.Get_Rendering as opsplt
import openseespy.postprocessing.ops_vis as opsv
from viktor.result import DownloadResult
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from pathlib import Path

def concrete_func(design_type, ecc_frp, fcc_retrofit, ax1, ay1, depth,width, cover, long_dia, trans_dia, nx, ny, n_leg_x, n_leg_y, coefficient,concrete_strength, steel_strength, s):
    
    # Inputs for concrete material analysis
    total_rebar = 2*nx + 2*(ny-2)
    d = depth-cover # clear height
    b_clear = width-cover # clear height
    A_long = 3.14*long_dia**2/4 #
    A_tra = 3.14*trans_dia**2/4 
    As = A_long*total_rebar
    bo = d-cover-trans_dia
    ho = b_clear-cover-trans_dia
    A = ho*bo/1000000
    A_tra_x = A_tra*n_leg_x
    A_tra_y = A_tra*n_leg_y
    ratio_x = (n_leg_x*bo*A_tra)/(ho*bo*s)
    ratio_y = (n_leg_y*ho*A_tra)/(ho*bo*s)
    total_ratio = ratio_x + ratio_y      

    # To get total k value help with distance between longitidunal rebars that 
    # supported by stirrups
    if ny == 2:
        total_ax = ax1**2
    elif ny == 3:
        total_ax = 2*(ax1**2)
    elif ny == 4:
        total_ax = 3*(ax1**2)
    elif ny == 5:
        total_ax = 4*(ax1**2)
    elif ny == 6:
        total_ax = 5*(ax1**2)
    elif ny == 7:
        total_ax = 6*(ax1**2)
    if nx == 2:
        total_ay = ay1**2
    elif nx == 3:
        total_ay = 2*(ay1**2)
    elif nx == 4:
        total_ay = 3*(ay1**2)
    elif nx == 5:
        total_ay = 4*(ay1**2)
    elif nx == 6:
        total_ay = 5*(ay1**2)
    elif nx == 7:
        total_ay = 6*(ay1**2)
    
    total_a = 2*total_ay+ 2*total_ax
    
    ke = (1-total_a/(6*bo*ho))*(1-s/(2*bo))*(1-s/(2*ho))*(1-As/(bo*ho))**-1

    if coefficient == "Nominal":
        fyd = steel_strength
    elif coefficient == "Expected":
        fyd = steel_strength*1.2
        
    fyh = steel_strength
    fex = ke*ratio_x*fyh
    fey = ke*ratio_y*fyh
    f1 = (fex+fey)/2
    if coefficient == "Nominal":
        fco = concrete_strength
    elif coefficient == "Expected":
        fco = concrete_strength*1.3
    fco = concrete_strength
    lambda_c = 2.254*sqrt(1+7.94*f1/fco)-2*(f1/fco)-1.254
    
    fcc_existing = float(format(fco*lambda_c, ".2f"))
    fsp = 0
    eco = 0.002
    ecu = 0.0035
    esp = 0.005
    
    if design_type == "Existing":
        ecc = eco*(1+5*(lambda_c-1))
        fcc = fcc_existing

    elif design_type == "Retrofit": 
        ecc = ecc_frp
        fcc = int(fcc_retrofit)

    
    # Ec = 5000*sqrt(fco)
    young_modulus_concrete = 5000*sqrt(concrete_strength)
    Esec = fcc/ecc
    Esec_unc = fco/eco
    r = young_modulus_concrete/(young_modulus_concrete-Esec)
    r_unc = young_modulus_concrete/(young_modulus_concrete-Esec_unc)
    f_cu = fco*(ecu/eco)*r_unc/(r_unc-1+(ecu/eco)**r_unc)
    
    e = 0
    fc_conf= []
    fc_unconf= []
    ec_conf = []
    ec_unconf = []
    x_conf = []
    x_unconf = []
    while e < 0.02:
            x = e/ecc
            fc = (fcc*x*r)/(r-1+x**r)
            fc_conf.append(format(fc, ".2f"))
            ec_conf.append(format(e, ".5f"))
            x_conf.append(format(x, ".3f"))
            e = e + 0.0001
            
    e = 0
    while e <= esp:
            x = e/eco
            if e <= ecu:
                fc = fco*x*r_unc/(r_unc-1+x**r_unc)
            elif e > ecu and e<=esp:
                fc = f_cu+(e-ecu)*((fsp-f_cu)/(esp-ecu))
            fc_unconf.append(format(fc, ".2f"))
            ec_unconf.append(format(e, ".5f"))
            x_unconf.append(format(x, ".3f"))
            e = e + 0.0001
            
    df_conf = pd.DataFrame(list(zip(ec_conf, fc_conf)), columns =['Strain', 'Stress'], dtype = float)
    
    we = min(ratio_x, ratio_y)*ke*(fyh/fco)
    
    if 0.0035+0.04*sqrt(we) <= 0.018:
        ecu_maximum = float(format(0.0035+0.04*sqrt(we), ".3f"))
    else:
        ecu_maximum = 0.018
        
    fcc_max = df_conf.loc[df_conf['Stress']==fcc]
    eco_max = fcc_max['Strain'].iloc[0]
    ecu_max = df_conf.loc[df_conf['Strain']==ecu_maximum]
    fcu_max = ecu_max['Stress'].iloc[0]
     
   
    return fcc, eco, esp, ecu_max, fcu_max, df_conf, eco_max, ecu_maximum, ecc

def MomentCurvature(name, concrete_output,steel_output,moment_curvature_output,secTag, b1, b2, h1, h2, axialLoad, maxK, numIncr):
    
    # Define two nodes at (0,0)
    node(1, 0.0, 0.0)
    node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    fix(1, 1, 1, 1)
    fix(2, 0, 1, 0)
    
    # Define element
    #                             tag ndI ndJ  secTag
    element('zeroLengthSection',  1,   1,   2,  secTag)      
         
    # Create recorder  
    recorder('Element', '-file', concrete_output, '-precision', int(5), '-time', '-dT', float(0.1) ,'-ele', 1, 'section', 'fiber', str(b1), str(h1), '1', 'stressStrain')
    recorder('Element', '-file', steel_output, '-precision', int(5), '-time', '-dT', float(0.1) ,'-ele', 1, 'section', 'fiber', str(b2), str(h2), '3', 'stressStrain')
    recorder('Node', '-file', moment_curvature_output, '-time','-dT', float(0.1) , '-node', 2, '-dof', 3, 'disp')

    # Define constant axial load
    timeSeries('Constant', 1)
    pattern('Plain', 1, 1)
    load(2, int(axialLoad), 0.0, 0.0)

    # Define analysis parameters
    integrator('LoadControl', 0.0)
    system('SparseGeneral', '-piv')
    test('NormUnbalance', 1e-6, 10)
    numberer('Plain')
    constraints('Plain')
    algorithm('Newton')
    analysis('Static')

    # Do one analysis for constant axial load
    analyze(1)

    # Define reference moment
    timeSeries('Linear', 2)
    pattern('Plain',2, 2)
    load(2, 0.0, 0.0, 1.0)

    # Compute curvature increment
    dK = maxK / numIncr

    # Use displacement control at node 2 for section analysis
    integrator('DisplacementControl', 2,3,dK,1,dK,dK)

    # Do the section analysis
    analyze(numIncr)
        
    
class LogoController(ViktorController):
    label = 'Logo design'
    parametrization = LogoParametrization

    viktor_convert_entity_field = True
   
    @staticmethod
    def interaction_diagram(params):
        
        #inputs
        analysis_type = params.tab_1.geometry.sectiontype
        concrete_strength = int(params.tab_1.material.concrete)
        steel_strength = int(params.tab_1.material.steel)
        gc = int(params.tab_1.material.gamma_concrete)
        gs = int(params.tab_1.material.gamma_steel)
        long_dia = params.tab_1.longitudinal.diameter
        trans_dia = params.tab_1.transverse.diameter
        b = params.tab_1.geometry.depth
        h = params.tab_1.geometry.width
        cover = params.tab_1.geometry.cover
        d = b-cover
        axial = params.tab_1.force.axialforce
        s = int(params.tab_1.transverse.space)
        nx = params.tab_1.longitudinal.nx
        ny = params.tab_1.longitudinal.ny
        length = params.tab_1.geometry.length
        mb = params.tab_1.shearforce.moment_bottom
        mt = params.tab_1.shearforce.moment_bottom
        ve = params.tab_1.shearforce.earthquake
        v_gqe = params.tab_1.shearforce.combination
        v_gqde = params.tab_1.shearforce.combination_d 
        degree = params.tab_1.geometry.degree
        sectiontype = params.tab_3.axial_frp.sectiontype
        rc = params.tab_3.axial_frp.rc 
        Ef = params.tab_3.axial_frp.youngmodulusofwrap
        ef = params.tab_3.axial_frp.strainofwrap
        eu = params.tab_3.axial_frp.ultimatestrainofwrap
        nf = params.tab_3.axial_frp.numberoffrp
        strip_type = params.tab_3.axial_frp.striptype
        tf = params.tab_3.axial_frp.thicknessofwrap
        sf = params.tab_3.axial_frp.distanceofwraps
        
        if strip_type == "Uncontinuous":
            fi_f = nf*wf*(b+h)*2*tf/(b*h*sf)
        else:
            fi_f = nf*(b+h)*2*tf/(b*h)
            
        if ef <= 0.004:
            ef_1 = ef
        else:
            ef_1 = 0.004
            
        ef_2 = 0.50*eu
        ef = min(ef_1,ef_2)
        
        if sectiontype == "Rectangular Section":
            kappa_a = 1-((b-2*rc)**2 + (h+2*rc)**2)/(3*b*h)
        elif sectiontype == "Circular Section":
            kappa_a = 1
        else:
            kappa_a = b/h
            
  
        f1 = 0.5*fi_f*kappa_a*ef*Ef
        fcc_1 = concrete_strength*(1+2.4*(f1/concrete_strength))
        fcc_2 = 1.2*concrete_strength
        fcc_frp = min(fcc_1,fcc_2) # Compression strength of the concrete to calculate the axial force
        
        concrete_list = [concrete_strength, fcc_frp]
        
        for each in concrete_list:
            #material
            young_modulus_concrete = 5000*sqrt(each)    
            fctk = 0.35*sqrt(each)
            fctd = fctk/gc
            fcd = each/gc
            fyd = steel_strength/gs
            young_modulus_steel = 200000
            s_confinement = max(min(min(b,h)/3,150,long_dia*6),50)
            a = trans_dia*25
            
            diameter_area = long_dia*long_dia*pi/4
            
            x_dir_list = []
            
            if degree == 0:
                d = h - cover
                k = 1
                x_reb = (b-2*cover)/(nx-1)
                increase = (b-2*cover)/(nx-1)
                while k <= nx:
                    if k == 1:
                        x_ = cover
                    elif k == nx:
                        x_ = b - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
            elif degree == 90:
                d = h - cover
                k = 1
                x_reb = (h-2*cover)/(ny-1)
                increase = (h-2*cover)/(ny-1)
                while k <= ny:
                    if k == 1:
                        x_ = cover
                    elif k == ny:
                        x_ = h - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
    
            ecu = 0.003
            esy = fyd/young_modulus_steel
            k1 = 0.82
            index_0 = 15
            index_90 = 15
            k =min(b,h)/index_0
            sensibility_0 = min(b,h)/index_0
            sensibility_90 = min(b,h)/index_90
     
            if degree == 0:
                cb = ecu*(h-cover)/(esy+ecu)
        
                if cb != 0:
                    if cb*k1 > h:
                        k1cb_ = h
                    else:
                        k1cb_ = cb*k1
                total_step = int(ceil((h-cover)/sensibility_0)+2)
                k1cb_list = []
                c_list = []
                deneme_mb_list = []
                c_index = 1
                total_layer_index = 1
                Force_List = []
                Moment_List = []
                while c_index <= total_step:
                    if c_index == 1:
                        c = h*1.3
                    elif c_index == 2:
                        c = h
                    elif c_index == total_step:
                        c = cover
                    else:
                        c = h - sensibility_0   
                        sensibility_0 = sensibility_0+k
                    if c*k1 > h:
                        k1cb = h 
                    else:
                        k1cb = c*k1 
                    k1cb_list.append(k1cb)
            
                    
                    rebar_number_list = []
                    rebar_distance_list = []
                    diameter_list = []
                    area_list = []
                    Le1_list = []
                    ec1_list = []
                    sigma_s_list = []
                    Fs1_list = []
                    Mb_list = []
                    distance_dif_cons = (b-2*cover)/(nx-1)
                    distance_dif = (b-2*cover)/(nx-1)
                    total_layer_index = 1
                    while total_layer_index <= nx:
                        if total_layer_index == 1:
                            rebar_number = ny 
                            rebar_distance = cover
                            area = rebar_number*diameter_area
                           
                        elif total_layer_index > 1 and total_layer_index < nx:
                            rebar_number = 2
                            rebar_distance = cover + distance_dif
                            area = rebar_number*diameter_area
                            distance_dif = distance_dif + distance_dif_cons
                            
                        else:
                            rebar_number = ny
                            rebar_distance = b-cover
                            area = rebar_number*diameter_area
                            
                        Le1 = c - rebar_distance    
                        if (Le1/c)*ecu > 0:
                            if (Le1/c)*ecu > esy:
                                es1 = esy
                            else:
                                es1 = (Le1/c)*ecu
                        elif (Le1/c)*ecu <= 0:
                            if (Le1/c)*ecu < esy*-1:
                                es1 = esy*-1
                            else:
                                es1 = (Le1/c)*ecu
                                
                        if es1*young_modulus_steel > 0:
                            if es1*young_modulus_steel > fyd:
                                sigma = fyd
                            else:
                                sigma = es1*young_modulus_steel
                        elif es1*young_modulus_steel <= 0:
                            if es1*young_modulus_steel <= fyd*-1:
                                sigma = fyd*-1
                            else:
                                sigma = es1*young_modulus_steel
                        Fc = 0.85*fcd*b*k1cb/1000
                              
                        Fs1 = (area*sigma/1000) 
                        Mb_ = Fc*(h/2-k1cb/2)/1000
                        Mb = Fs1*(h/2-rebar_distance)/1000  
                        
                        rebar_number_list.append(rebar_number)
                        rebar_distance_list.append(rebar_distance)
                        diameter_list.append(long_dia)
                        area_list.append(area)
                        Le1_list.append(int(Le1))
                        ec1_list.append(es1)
                        sigma_s_list.append(sigma)
                        Fs1_list.append(Fs1)
                        Mb_list.append(Mb)
                        if c_index == 1:
                            deneme_mb_list.append(Mb)
                
                 
                        #print(Fs1_list)
                        total_layer_index = total_layer_index +1
                    total_Fc = int(ceil(sum(Fs1_list) + Fc))
                    total_Mb = int(ceil(sum(Mb_list) + Mb_))
                    total_area = int(ceil(sum(area_list)))
                    Moment_List.append(total_Mb)
                    Force_List.append(total_Fc)
                    c_list.append(c)
                    c_index = c_index +1 
                    
            elif degree == 90:
                cb = ecu*(b-cover)/(esy+ecu)
                
                if cb != 0:
                    if cb*k1 > b:
                        k1cb_ = b
                    else:
                        k1cb_ = cb*k1
                total_step = int(ceil((b-cover)/sensibility_0)+2)
                k1cb_list = []
                c_list = []
                deneme_mb_list = []
                c_index = 1
                total_layer_index = 1
                Force_List = []
                Moment_List = []
                while c_index <= total_step:
                    if c_index == 1:
                        c = b*1.3
                    elif c_index == 2:
                        c = b
                    elif c_index == total_step:
                        c = cover
                    else:
                        c = b - sensibility_0   
                        sensibility_0 = sensibility_0+k
                    if c*k1 > b:
                        k1cb = b 
                    else:
                        k1cb = c*k1 
                    k1cb_list.append(k1cb)
            
                    
                    rebar_number_list = []
                    rebar_distance_list = []
                    diameter_list = []
                    area_list = []
                    Le1_list = []
                    ec1_list = []
                    sigma_s_list = []
                    Fs1_list = []
                    Mb_list = []
                    distance_dif_cons = (h-2*cover)/(ny-1)
                    distance_dif = (h-2*cover)/(ny-1)
                    total_layer_index = 1
                    while total_layer_index <= ny:
                        if total_layer_index == 1:
                            rebar_number = nx 
                            rebar_distance = cover
                            area = rebar_number*diameter_area
                           
                        elif total_layer_index > 1 and total_layer_index < ny:
                            rebar_number = 2
                            rebar_distance = cover + distance_dif
                            area = rebar_number*diameter_area
                            distance_dif = distance_dif + distance_dif_cons
                            
                        else:
                            rebar_number = nx
                            rebar_distance = h-cover
                            area = rebar_number*diameter_area
                            
                        Le1 = c - rebar_distance    
                        if (Le1/c)*ecu > 0:
                            if (Le1/c)*ecu > esy:
                                es1 = esy
                            else:
                                es1 = (Le1/c)*ecu
                        elif (Le1/c)*ecu <= 0:
                            if (Le1/c)*ecu < esy*-1:
                                es1 = esy*-1
                            else:
                                es1 = (Le1/c)*ecu
                                
                        if es1*young_modulus_steel > 0:
                            if es1*young_modulus_steel > fyd:
                                sigma = fyd
                            else:
                                sigma = es1*young_modulus_steel
                        elif es1*young_modulus_steel <= 0:
                            if es1*young_modulus_steel <= fyd*-1:
                                sigma = fyd*-1
                            else:
                                sigma = es1*young_modulus_steel
                        Fc = 0.85*fcd*h*k1cb/1000
                              
                        Fs1 = (area*sigma/1000) 
                        Mb_ = Fc*(b/2-k1cb/2)/1000
                        Mb = Fs1*(b/2-rebar_distance)/1000  
                        
                        rebar_number_list.append(rebar_number)
                        rebar_distance_list.append(rebar_distance)
                        diameter_list.append(long_dia)
                        area_list.append(area)
                        Le1_list.append(int(Le1))
                        ec1_list.append(es1)
                        sigma_s_list.append(sigma)
                        Fs1_list.append(Fs1)
                        Mb_list.append(Mb)
                        if c_index == 1:
                            deneme_mb_list.append(Mb)
                
                 
                        #print(Fs1_list)
                        total_layer_index = total_layer_index +1
                    total_Fc = int(ceil(sum(Fs1_list) + Fc))
                    total_Mb = int(ceil(sum(Mb_list) + Mb_))
                    total_area = int(ceil(sum(area_list)))
                    Moment_List.append(total_Mb)
                    Force_List.append(total_Fc)
                    c_list.append(c)
                    c_index = c_index +1 
            
            Nc = (b*h*0.85*fcd+fyd*total_area)/1000
            Nt = -1*fyd*total_area/1000
            Moment_List.insert(0, 0)
            Moment_List.insert(total_step+1, 0)
            Force_List.insert(0, Nc)
            Force_List.insert(total_step+1, Nt)
            
            if each == concrete_strength:
                moment_existing = Moment_List
                force_existing = Force_List
            else:
                moment_retrofit = Moment_List
                force_retrofit= Force_List
                
        return moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp
    
    @PlotlyView("Interaction Diagram", duration_guess=1)
    def get_data_view1(self, params, **kwargs):

        analysis_type = params.tab_1.geometry.sectiontype
        moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp = self.interaction_diagram(params)
        
        df_m_n = pd.DataFrame(list(zip(force_existing, moment_existing)), columns =['Axial Force', 'Moment'], dtype = float)
        moment_existing_180 = [x * -1 for x in moment_existing]
        df_m_n_180 = pd.DataFrame(list(zip(force_existing, moment_existing_180)), columns =['Axial Force', 'Moment'], dtype = float)
        moment_retrofit_180 = [x * -1 for x in moment_retrofit]
    
        if analysis_type == "Retrofit":
            fig = go.Figure()
    
            fig.add_trace(go.Scatter(x=moment_existing, y=force_existing,  
                                mode='lines+markers',
                                name='0-Degree_Existing'))
            fig.add_trace(go.Scatter(x=moment_existing_180, y=force_existing,
                                mode='lines+markers',
                                name='180-Degree_Existing'))
            fig.add_trace(go.Scatter(x=moment_retrofit, y=force_retrofit,
                                mode='lines+markers',
                                name='180-Degree_Retrofit'))
            fig.add_trace(go.Scatter(x=moment_retrofit_180, y=force_retrofit,
                                mode='lines+markers',
                                name='180-Degree_Retrofit'))
        else:
            fig = go.Figure()
    
            fig.add_trace(go.Scatter(x=moment_existing, y=force_existing,  
                                mode='lines+markers',
                                name='0-Degree_Existing'))
            fig.add_trace(go.Scatter(x=moment_existing_180, y=force_existing,
                                mode='lines+markers',
                                name='180-Degree_Existing'))
        
        return PlotlyResult(fig.to_json())
    
    @staticmethod
    def calculate_moment_and_curvature(params):
        
        moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp = LogoController.interaction_diagram(params)
        concrete_strength = int(params.tab_1.material.concrete)
        steel_strength = int(params.tab_1.material.steel)
        gc = int(params.tab_1.material.gamma_concrete)
        gs = int(params.tab_1.material.gamma_steel)
        long_dia = params.tab_1.longitudinal.diameter
        trans_dia = params.tab_1.transverse.diameter
        b = params.tab_1.geometry.depth
        h = params.tab_1.geometry.width
        depth = params.tab_1.geometry.depth
        width = params.tab_1.geometry.width
        cover = params.tab_1.geometry.cover
        d = b-cover
        axialLoad = params.tab_1.force.axialforce
        s = int(params.tab_1.transverse.space)
        nx = params.tab_1.longitudinal.nx
        ny = params.tab_1.longitudinal.ny
        length = params.tab_1.geometry.length
        mb = params.tab_1.shearforce.moment_bottom
        mt = params.tab_1.shearforce.moment_bottom
        ve = params.tab_1.shearforce.earthquake
        v_gqe = params.tab_1.shearforce.combination
        v_gqde = params.tab_1.shearforce.combination_d 
        degree = params.tab_1.geometry.degree
        n_leg_x = params.tab_1.transverse.nlegx
        n_leg_y = params.tab_1.transverse.nlegy
        ax1 = params.tab_2.concreteproperties.ax
        ay1 = params.tab_2.concreteproperties.ay
        coefficient = str(params.tab_1.material.coefficient)
        s = params.tab_1.transverse.space
        eu = params.tab_3.axial_frp.ultimatestrainofwrap
        analysis_type = params.tab_1.geometry.sectiontype
        confinement_type = params.tab_1.transverse.confinementtype
        
        #material
        young_modulus_concrete = 5000*sqrt(concrete_strength)    
        fctk = 0.35*sqrt(concrete_strength)
        fctd = fctk/gc
        fcd = concrete_strength/gc
        fyd = steel_strength/gs
        young_modulus_steel = 200000
        s_confinement = max(min(min(b,h)/3,150,long_dia*6),50)
        a = trans_dia*25
        
        ef_1 = 0.50*eu            
        ef_2 = 0.01
        ef = min(ef_1,ef_2)
        
        ecc_frp = 0.002*(1+15*((f1/fcd)**0.75))
        fcc_retrofit = fcc_frp
        
        diameter_area = long_dia*long_dia*pi/4

        design_list = ["Existing", "Retrofit"]
        
        for each in design_list:
            design_type = each
            fcc, eco, esp, ecu_max, fcu_max, df_conf, eco_max, ecu_maximum, ecc = concrete_func(design_type, ecc_frp, fcc_retrofit, ax1, ay1, depth,width, cover, long_dia, trans_dia, nx, ny, n_leg_x, n_leg_y, coefficient, concrete_strength, steel_strength, s)
            x_dir_list = []
    
            if degree == 0:
                d = h - cover
                k = 1
                x_reb = (b-2*cover)/(nx-1)
                increase = (b-2*cover)/(nx-1)
                while k <= nx:
                    if k == 1:
                        x_ = cover
                    elif k == nx:
                        x_ = b - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
            elif degree == 90:
                d = h - cover
                k = 1
                x_reb = (h-2*cover)/(ny-1)
                increase = (h-2*cover)/(ny-1)
                while k <= ny:
                    if k == 1:
                        x_ = cover
                    elif k == ny:
                        x_ = h - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
            
            ecu = 0.002
            esy = fyd/young_modulus_steel
            
            if degree == 0:
                y1 = int(b/2.0)
                z1 = int(h/2.0)
            
                length_le1 = len(x_dir_list)          
                
            elif degree == 90:
                y1 = int(h/2.0)
                z1 = int(b/2.0)
                length_le1 = len(x_dir_list)          
    
            wipe()
            print("Start MomentCurvature.py example")
            
            concrete_output = "concrete_strain.txt"
            steel_output = "steel_strain.txt"
            moment_curvature_output = "moment_curvature.txt"   
            name = "name"
            
            # ------------------------------ Moment Curvature ------------------------------
            # Define model builder
            # --------------------
            model('basic','-ndm',2,'-ndf',3)                   
            secTag = 1        
            
            # Cover concrete (unconfined)
            # uniaxialMaterial('Concrete04',2, -fcd,  -eco,  -esp,  young_modulus_concrete, 0.0, 0.0, 0,1)
            # uniaxialMaterial('Concrete01',1, -fcc, -eco_max,  -fcu_max,  -ecu_maximum)
            uniaxialMaterial('Concrete02',1, -fcu_max, -ecu,  -0.2*fcu_max,  -0.02, 0.1, 0.1*fcu_max, (0.1*fcu_max)/ecu)
            # uniaxialMaterial('Concrete01',1, -fcc, -eco_max,  -fcu_max,  -0.02)
            # uniaxialMaterial('Concrete04',1, float(fcc),  ecc,  -0.02,  Ec)
            
            # # Cover concrete (unconfined)
            # uniaxialMaterial('Concrete04',2, float(fc0),  -0.002,  -0.004,  Ec)
            uniaxialMaterial('Concrete02',2, -fcd, -ecu,  -0.2*fcd,  -0.02, 0.1, 0.1*fcd, (0.1*fcd)/ecu)
            # uniaxialMaterial('Concrete01',2, -fcd,  -ecu,  0.0,  -esp)
            
            # STEEL
            # Reinforcing steel 
            by = 0.01
            
            #                        tag  fy E0    b
            uniaxialMaterial('Steel01', 3, int(fyd), young_modulus_steel, by)
            
            # Define cross-section for nonlinear columns
            # ------------------------------------------
            
            if degree == 0:
                width = h
                depth = b
                cover = cover
                number_of_layer = length_le1
                n_top = ny
                n_bot = ny
                n_int = 2
                
                
                b1 = width/2 - cover
                b2 = (width/2 - cover)*-1
                h1 = depth/2 - cover
                h2 = (depth/2 - cover)*-1
                
                # some variables derived from the parameters
                y1 = depth/2.0
                z1 = width/2.0
                total_y = depth - 2*cover
                total_y_layer = total_y/(number_of_layer-1)
                total_y_layer_step = total_y/(number_of_layer-1)
                
                section('Fiber', 1)
                
                # Create the concrete core fibers
                if confinement_type == "Confined" or confinement_type == "Unconfined" and design_type == "Retrofit":
                    patch('rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                    
                    
                    # Create the concrete cover fibers (top, bottom, left, right)
                    patch('rect',2,20,1 ,-y1, z1-cover, y1, z1)
                    patch('rect',2,20,1 ,-y1, -z1, y1, cover-z1)
                    patch('rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover)
                    patch('rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover)

                
                    top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                    bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                    
                    fib_sec_1 = [['section', 'Fiber', 1],
                    ['patch', 'rect',2,20,1 ,-y1, z1-cover, y1, z1],
                    ['patch', 'rect',2,20,1 ,-y1, -z1, y1, cover-z1],
                    ['patch', 'rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover],
                    ['patch', 'rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover],
                    ['patch', 'rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                    top,
                    bottom]
                    
                elif confinement_type == "Unconfined" and design_type == "Existing":
                    patch('rect',2,10,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                
                    top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                    bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                    
                    fib_sec_1 = [['section', 'Fiber', 1],
                    ['patch', 'rect',2,10,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                    top,
                    bottom]
                
                elif confinement_type == "Confined" and design_type == "Existing":
                    patch('rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                    
                    
                    # Create the concrete cover fibers (top, bottom, left, right)
                    patch('rect',2,20,1 ,-y1, z1-cover, y1, z1)
                    patch('rect',2,20,1 ,-y1, -z1, y1, cover-z1)
                    patch('rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover)
                    patch('rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover)

                
                    top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                    bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                    
                    fib_sec_1 = [['section', 'Fiber', 1],
                    ['patch', 'rect',2,20,1 ,-y1, z1-cover, y1, z1],
                    ['patch', 'rect',2,20,1 ,-y1, -z1, y1, cover-z1],
                    ['patch', 'rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover],
                    ['patch', 'rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover],
                    ['patch', 'rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                    top,
                    bottom]
                    
                # plt.savefig('fibsec_rc.png')
                
                # # Create the reinforcing fibers (left, middle, right)
                
                layer('straight', 3, n_top, diameter_area, y1-cover, cover-z1, y1-cover, z1-cover)
                layer('straight', 3, n_bot, diameter_area, cover-y1, cover-z1, cover-y1, z1-cover)
                
                total_int_layer = number_of_layer-2
                int_layer = 1
                while int_layer <= total_int_layer:
                
                    layer('straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia))
                    int_layer_def = ['layer','straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia)]
                    fib_sec_1.append(int_layer_def)
                    total_y_layer = total_y_layer + total_y_layer_step
                    int_layer = int_layer +1
                    
                # d -- from cover to rebar
                d = depth-cover
            
            elif degree == 90:
                width = b
                depth = h
                cover = cover
                number_of_layer = length_le1
                n_top = nx
                n_bot = nx
                n_int = 2
                
                
                b1 = width/2 - cover
                b2 = (width/2 - cover)*-1
                h1 = depth/2 - cover
                h2 = (depth/2 - cover)*-1
                
                # some variables derived from the parameters
                y1 = depth/2.0
                z1 = width/2.0
                total_y = depth - 2*cover
                total_y_layer = total_y/(number_of_layer-1)
                total_y_layer_step = total_y/(number_of_layer-1)
                
                section('Fiber', 1)
                
                # Create the concrete core fibers
                patch('rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                
                
                # Create the concrete cover fibers (top, bottom, left, right)
                patch('rect',2,20,1 ,-y1, z1-cover, y1, z1)
                patch('rect',2,20,1 ,-y1, -z1, y1, cover-z1)
                patch('rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover)
                patch('rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover)
                
                top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                
                fib_sec_1 = [['section', 'Fiber', 1],
                ['patch', 'rect',2,20,1 ,-y1, z1-cover, y1, z1],
                ['patch', 'rect',2,20,1 ,-y1, -z1, y1, cover-z1],
                ['patch', 'rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover],
                ['patch', 'rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover],
                ['patch', 'rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                top,
                bottom]
                    
                
                # plt.savefig('fibsec_rc.png')
                
                # # Create the reinforcing fibers (left, middle, right)
                
                layer('straight', 3, n_top, diameter_area, y1-cover, cover-z1, y1-cover, z1-cover)
                layer('straight', 3, n_bot, diameter_area, cover-y1, cover-z1, cover-y1, z1-cover)
                
                total_int_layer = number_of_layer-2
                int_layer = 1
                while int_layer <= total_int_layer:
                
                    layer('straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia))
                    int_layer_def = ['layer','straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia)]
                    fib_sec_1.append(int_layer_def)
                    total_y_layer = total_y_layer + total_y_layer_step
                    int_layer = int_layer +1
                    
                # d -- from cover to rebar
                d = width-cover
                
            matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
            m_c = opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor) #Figure of the section
            
            plt.axis('equal')    
            
            # Estimate yield curvature
            # (Assuming no axial load and only top and bottom steel)
            # d -- from cover to rebar
            
            # steel yield strain
            epsy = fyd/young_modulus_steel
            Ky = epsy/(0.7*d)
            
            # Print estimate to standard output
            print("Estimated yield curvature: ", Ky)
               
            # Target ductility for analysis
            mu = 30
            
            # Number of analysis increments
            numIncr = 500
            
            # Call the section analysis procedure
            
            results = open('results.out','a+')
            
            # u = nodeDisp(2,3)
            # if abs(u-0.00190476190476190541)<1e-12:
            #     results.write('PASSED : MomentCurvature.py\n');
            #     print("Passed!")
            # else:
            #     results.write('FAILED : MomentCurvature.py\n');
            #     print("Failed!")
            
            results.close()
            
            print("==========================")
            
            axialDesign = -1000*axialLoad
            
            MomentCurvature(name, concrete_output,steel_output,moment_curvature_output, secTag, b1, b2, h1, h2, axialDesign, Ky*mu, numIncr)
            
            
            # Reading of Moment Curvature Results
            with open(moment_curvature_output) as f:
                coords = f.read().split()
            
            # Splitting data as Moment & Curvature
            moment = coords[0::2]
            curvature = coords[1::2]
            
            moment = moment[:-1]
            curvature = curvature[:-1]
            
            df = pd.DataFrame(list(zip(curvature, moment)), columns =['Curvature', 'Moment'], dtype = float)
            df['Curvature'] = 1000*df['Curvature']
            df['Moment'] = df['Moment']/1000000
            df.plot(kind='line',x='Curvature',y='Moment',color='red')
            moment_list = df['Moment'].tolist()
            curvature_list = df["Curvature"].tolist()
            
            if each == "Existing":
                moment_existing = moment_list
                curvature_existing = curvature_list
                df_existing = pd.DataFrame(list(zip(curvature_existing, moment_existing)), columns =['Curvature', 'Moment'], dtype = float)
            elif each == "Retrofit":
                moment_retrofit = moment_list
                curvature_retrofit= curvature_list
                df_retrofit = pd.DataFrame(list(zip(curvature_retrofit, moment_retrofit)), columns =['Curvature', 'Moment'], dtype = float)
            
            
        
        def convert_design_m_c(df):
            return df.to_csv().encode('utf-8')

        csv_design = convert_design_m_c(df)
        document_code_mc = "Moment - Curvature (0 Degree) - Press to Download"
        csv_file_mc = "m_n_interaction_0.csv"
        
            
        return moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing
    
    @PlotlyView("Moment Curvature", duration_guess=1)
    def get_plotly_view(self, params, **kwargs):
        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing = self.calculate_moment_and_curvature(params)
        
        fig = go.Figure()
        # fig = px.line(x=Moment_List, y=Force_List, labels={'x':'Moment', 'y':'Axial Force'})
        # Add traces
        # fig.add_trace(go.Scatter(x=moment_list, y=curvature_list,  
        #                     mode='lines+markers',
        #                     name='0-Degree'))
        # # fig.add_trace(go.Scatter(x=curvature_list, y=moment_list,
        # #                     mode='lines+markers',
        # #                     name='180-Degree'))
        

        fig = go.Figure()

        fig.add_trace(go.Scatter(x=curvature_existing, y=moment_existing,  
                            mode='lines+markers',
                            name='0-Degree_Existing'))
        fig.add_trace(go.Scatter(x=curvature_retrofit, y=moment_retrofit,
                            mode='lines+markers',
                            name='0-Degree_Retrofit'))
        
        #fig = px.line(x=curvature_list, y=moment_list, labels={'x':'Curvature', 'y':'Moment'})
        
        return PlotlyResult(fig.to_json())
    

    @GeometryView('3D', duration_guess=1)
    def get_3d_view(self, params, **kwargs):

        depth = params.tab_1.geometry.depth
        width = params.tab_1.geometry.width
        height = params.tab_1.geometry.length

        beam = SquareBeam(depth, width, height)

        # axis_system_location = Point(2, 2, 2)
        # axis_system = CartesianAxes(axis_system_location, axis_length=0.5, axis_diameter=0.5)

        geometries = [beam]

        return GeometryResult(geometries)
    
    def download_mass_breakdown_existing(self, params, **kwargs):

        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing = self.calculate_moment_and_curvature(params)

        # df = pd.DataFrame(moment_list)
        
        return DownloadResult(df_existing.to_csv(), 'moment_curvature_existing.csv')
    
    def download_mass_breakdown_retrofit(self, params, **kwargs):

        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing = self.calculate_moment_and_curvature(params)

        # df = pd.DataFrame(moment_list)
        
        return DownloadResult(df_retrofit.to_csv(), 'moment_curvature_retrofit.csv')
    
    @DataView('Shear Capacity', duration_guess=1)
    def get_data_view(self, params, **kwargs):

        #inputs
        moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp = self.interaction_diagram(params)
        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum = self.calculate_moment_and_curvature(params)
        concrete_strength = int(params.tab_1.material.concrete)
        steel_strength = int(params.tab_1.material.steel)
        gc = int(params.tab_1.material.gamma_concrete)
        gs = int(params.tab_1.material.gamma_steel)
        long_dia = params.tab_1.longitudinal.diameter
        trans_dia = params.tab_1.transverse.diameter
        depth = params.tab_1.geometry.depth
        width = params.tab_1.geometry.width
        cover = params.tab_1.geometry.cover
        d = depth-cover
        axial = params.tab_1.force.axialforce
        s = int(params.tab_1.transverse.space)
        n_leg_x = params.tab_1.transverse.nlegx
        n_leg_y = params.tab_1.transverse.nlegy
        length = params.tab_1.geometry.length
        mb = params.tab_1.shearforce.moment_bottom
        mt = params.tab_1.shearforce.moment_bottom
        ve = params.tab_1.shearforce.earthquake
        v_gqe = params.tab_1.shearforce.combination
        v_gqde = params.tab_1.shearforce.combination_d 
        degree = params.tab_1.geometry.degree
        nf = params.tab_3.frp.numberoffrp
        tf = params.tab_3.frp.thicknessofwrap
        wf = params.tab_3.frp.widthofwrap
        Ef = params.tab_3.frp.youngmodulusofwrap
        ef = params.tab_3.frp.strainofwrap
        d_wrap = params.tab_3.frp.effectivedepth
        sf = params.tab_3.frp.distanceofwraps
        strip_type = params.tab_3.frp.striptype
        eu = params.tab_3.frp.ultimatestrainofwrap
    
        #material
        young_modulus_concrete = 5000*sqrt(concrete_strength)    
        fctk = 0.35*sqrt(concrete_strength)
        fctd = fctk/gc
        fcd = concrete_strength/gc
        fyd = steel_strength/gs
        young_modulus_steel = 200000
        s_confinement = max(min(min(depth,width)/3,150,long_dia*6),50)
        a = trans_dia*25
        
        if axial < 0:
            gamma = -0.3
        else:
            gamma = 0.07

        #geometry
        
        axial_limit = depth*width*0.2*concrete_strength/1000
        axial_limit_ts500 = depth*width*0.05*concrete_strength/1000
        depth_clear = depth-2*cover-trans_dia
        width_clear = width-2*cover-trans_dia
        area_clear = depth_clear*width_clear
        area = depth*width
        trans_area_x = pi*trans_dia*trans_dia*n_leg_x/4 
        trans_area_y = pi*trans_dia*trans_dia*n_leg_y/4 
        s_mid = min(min(depth,width)/2,200)
        
        #shear force
        ve_capacity = (mb+mt)/length
        
        ve_final = min(ve_capacity, v_gqde)
        
        if axial < axial_limit_ts500 and ve > v_gqe/2:
            vc = 0
        else:
            vc = ((0.65*fctd*d*depth)*(1+gamma*axial*1000/(depth*d))*0.001)*0.8
        vw = trans_area_x*fyd*d/(s*1000)
        vr = vc + vw
        
        v_tbdy2018 = 0.85*depth*d*sqrt(fcd)*0.001
        
        if vr > ve_final and v_tbdy2018 > ve_final:
            shear_check = "TS500 and TBDY 2018 check are ok!"
        elif vr > ve_final and v_tbdy2018 < ve_final:
            shear_check = "TS500 check is ok - TBDY 2018 check is not ok!"
        elif vr < ve_final and v_tbdy2018 > ve_final:
            shear_check = "TS500 check is not ok - TBDY 2018 check is ok!"   
        elif vr < ve_final and v_tbdy2018 < ve_final:
            shear_check = "TS500 and TBDY 2018 check are not ok!" 

        
        if axial > axial_limit:
            con_reb_area_x = max(0.30*s_confinement*depth_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*depth_clear*(fcd/fyd))
            con_reb_area_y = max(0.30*s_confinement*width_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*width_clear*(fcd/fyd))
        else:
            con_reb_area_x = max(0.30*s_confinement*depth_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*depth_clear*(fcd/fyd))*2/3
            con_reb_area_y = max(0.30*s_confinement*width_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*width_clear*(fcd/fyd))*2/3
           
        if trans_area_x < con_reb_area_x:
            warning = "Transverse spacing is not enough for confinement zone!"
        else:
            warning = "Transverse spacing is enough for confinement zone."
        
        if ef <= 0.004:
            ef_1 = ef
        else:
            ef_1 = 0.004
            
        ef_2 = 0.50*eu
        ef = min(ef_1,ef_2)
         
        if strip_type == "Uncontinuous":
            sf = min(wf+d/4, sf)
            v_wrap = 2*nf*tf*wf*Ef*ef*d_wrap/sf/1000
        else:
            v_wrap = 2*nf*tf*Ef*ef*d_wrap/1000
        
        data = DataGroup(
            spacing_conf = DataItem(warning, s, suffix='mm', number_of_decimals=2),
            spacing_mid = DataItem("Transverse space is: ", s_mid, suffix='mm', number_of_decimals=2),
            tbdy_2018 = DataItem("TBDY2018 - Capacity:", v_tbdy2018, suffix='kN',number_of_decimals=2),
            # ef_final = DataItem("Unit Strain of FRP:", ef, number_of_decimals=3),
            # ecc_final = DataItem("Unit Strain of FRP:", ecc, number_of_decimals=3),
            # ecu_maximum1 = DataItem("Ecu Max:", ecu_maximum, number_of_decimals=3),
            # f1 = DataItem("f1:", f1, number_of_decimals=3),
            # fcc_frp = DataItem("fcc_frp:", fcc_frp, number_of_decimals=3),
            shearcapacity_stirrups = DataItem("Shear Force Capacity (Vw): ", vw, suffix='kN', number_of_decimals=2),
            shearcapacity_concrete = DataItem("Shear Force Capacity (Vc): ", vc, suffix='kN', number_of_decimals=2),
            frp_capacity = DataItem("Shear Capacity (Vf): ", v_wrap, suffix='kN', number_of_decimals=2),
            shearforce_capacity = DataItem("Total Shear Force Capacity without FRP (Vr): ", vr, suffix='kN', number_of_decimals=2),
            shearforce_capacity_frp = DataItem("Total Shear Force Capacity with FRP (Vf): ", v_wrap+vr, suffix='kN', number_of_decimals=2),
            
            # capacity_check = DataItem("Shear Capacity Check:", shear_check)

        )
        
        return DataResult(data)
    from math import cos, sin, pi, sqrt, ceil
import pandas as pd
from viktor.core import ViktorController
from viktor.geometry import Line, Point, RectangularExtrusion, Material, Point, CartesianAxes, SquareBeam
from viktor.views import GeometryResult, GeometryView, DataResult, DataGroup, DataView,DataView, DataItem, PlotlyView, PlotlyResult, PlotlyAndDataView, PNGResult  
from .parametrization import LogoParametrization, DownloadButton
import matplotlib.pyplot as plt
from io import StringIO
from viktor.views import SVGView, SVGResult, PNGView
from openseespy.opensees import *
import openseespy.postprocessing.Get_Rendering as opsplt
import openseespy.postprocessing.ops_vis as opsv
from viktor.result import DownloadResult
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from pathlib import Path

def concrete_func(design_type, ecc_frp, fcc_retrofit, ax1, ay1, depth,width, cover, long_dia, trans_dia, nx, ny, n_leg_x, n_leg_y, coefficient,concrete_strength, steel_strength, s):
    
    # Inputs for concrete material analysis
    total_rebar = 2*nx + 2*(ny-2)
    d = depth-cover # clear height
    b_clear = width-cover # clear height
    A_long = 3.14*long_dia**2/4 #
    A_tra = 3.14*trans_dia**2/4 
    As = A_long*total_rebar
    bo = d-cover-trans_dia
    ho = b_clear-cover-trans_dia
    A = ho*bo/1000000
    A_tra_x = A_tra*n_leg_x
    A_tra_y = A_tra*n_leg_y
    ratio_x = (n_leg_x*bo*A_tra)/(ho*bo*s)
    ratio_y = (n_leg_y*ho*A_tra)/(ho*bo*s)
    total_ratio = ratio_x + ratio_y      

    # To get total k value help with distance between longitidunal rebars that 
    # supported by stirrups
    if ny == 2:
        total_ax = ax1**2
    elif ny == 3:
        total_ax = 2*(ax1**2)
    elif ny == 4:
        total_ax = 3*(ax1**2)
    elif ny == 5:
        total_ax = 4*(ax1**2)
    elif ny == 6:
        total_ax = 5*(ax1**2)
    elif ny == 7:
        total_ax = 6*(ax1**2)
    if nx == 2:
        total_ay = ay1**2
    elif nx == 3:
        total_ay = 2*(ay1**2)
    elif nx == 4:
        total_ay = 3*(ay1**2)
    elif nx == 5:
        total_ay = 4*(ay1**2)
    elif nx == 6:
        total_ay = 5*(ay1**2)
    elif nx == 7:
        total_ay = 6*(ay1**2)
    
    total_a = 2*total_ay+ 2*total_ax
    
    ke = (1-total_a/(6*bo*ho))*(1-s/(2*bo))*(1-s/(2*ho))*(1-As/(bo*ho))**-1

    if coefficient == "Nominal":
        fyd = steel_strength
    elif coefficient == "Expected":
        fyd = steel_strength*1.2
        
    fyh = steel_strength
    fex = ke*ratio_x*fyh
    fey = ke*ratio_y*fyh
    f1 = (fex+fey)/2
    if coefficient == "Nominal":
        fco = concrete_strength
    elif coefficient == "Expected":
        fco = concrete_strength*1.3
    fco = concrete_strength
    lambda_c = 2.254*sqrt(1+7.94*f1/fco)-2*(f1/fco)-1.254
    
    fcc_existing = float(format(fco*lambda_c, ".2f"))
    fsp = 0
    eco = 0.002
    ecu = 0.0035
    esp = 0.005
    
    if design_type == "Existing":
        ecc = eco*(1+5*(lambda_c-1))
        fcc = fcc_existing

    elif design_type == "Retrofit": 
        ecc = ecc_frp
        fcc = int(fcc_retrofit)

    
    # Ec = 5000*sqrt(fco)
    young_modulus_concrete = 5000*sqrt(concrete_strength)
    Esec = fcc/ecc
    Esec_unc = fco/eco
    r = young_modulus_concrete/(young_modulus_concrete-Esec)
    r_unc = young_modulus_concrete/(young_modulus_concrete-Esec_unc)
    f_cu = fco*(ecu/eco)*r_unc/(r_unc-1+(ecu/eco)**r_unc)
    
    e = 0
    fc_conf= []
    fc_unconf= []
    ec_conf = []
    ec_unconf = []
    x_conf = []
    x_unconf = []
    while e < 0.02:
            x = e/ecc
            fc = (fcc*x*r)/(r-1+x**r)
            fc_conf.append(format(fc, ".2f"))
            ec_conf.append(format(e, ".5f"))
            x_conf.append(format(x, ".3f"))
            e = e + 0.0001
            
    e = 0
    while e <= esp:
            x = e/eco
            if e <= ecu:
                fc = fco*x*r_unc/(r_unc-1+x**r_unc)
            elif e > ecu and e<=esp:
                fc = f_cu+(e-ecu)*((fsp-f_cu)/(esp-ecu))
            fc_unconf.append(format(fc, ".2f"))
            ec_unconf.append(format(e, ".5f"))
            x_unconf.append(format(x, ".3f"))
            e = e + 0.0001
            
    df_conf = pd.DataFrame(list(zip(ec_conf, fc_conf)), columns =['Strain', 'Stress'], dtype = float)
    
    we = min(ratio_x, ratio_y)*ke*(fyh/fco)
    
    if 0.0035+0.04*sqrt(we) <= 0.018:
        ecu_maximum = float(format(0.0035+0.04*sqrt(we), ".3f"))
    else:
        ecu_maximum = 0.018
        
    fcc_max = df_conf.loc[df_conf['Stress']==fcc]
    eco_max = fcc_max['Strain'].iloc[0]
    ecu_max = df_conf.loc[df_conf['Strain']==ecu_maximum]
    fcu_max = ecu_max['Stress'].iloc[0]
     
   
    return fcc, eco, esp, ecu_max, fcu_max, df_conf, eco_max, ecu_maximum, ecc

def MomentCurvature(name, concrete_output,steel_output,moment_curvature_output,secTag, b1, b2, h1, h2, axialLoad, maxK, numIncr):
    
    # Define two nodes at (0,0)
    node(1, 0.0, 0.0)
    node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    fix(1, 1, 1, 1)
    fix(2, 0, 1, 0)
    
    # Define element
    #                             tag ndI ndJ  secTag
    element('zeroLengthSection',  1,   1,   2,  secTag)      
         
    # Create recorder  
    recorder('Element', '-file', concrete_output, '-precision', int(5), '-time', '-dT', float(0.1) ,'-ele', 1, 'section', 'fiber', str(b1), str(h1), '1', 'stressStrain')
    recorder('Element', '-file', steel_output, '-precision', int(5), '-time', '-dT', float(0.1) ,'-ele', 1, 'section', 'fiber', str(b2), str(h2), '3', 'stressStrain')
    recorder('Node', '-file', moment_curvature_output, '-time','-dT', float(0.1) , '-node', 2, '-dof', 3, 'disp')

    # Define constant axial load
    timeSeries('Constant', 1)
    pattern('Plain', 1, 1)
    load(2, int(axialLoad), 0.0, 0.0)

    # Define analysis parameters
    integrator('LoadControl', 0.0)
    system('SparseGeneral', '-piv')
    test('NormUnbalance', 1e-6, 10)
    numberer('Plain')
    constraints('Plain')
    algorithm('Newton')
    analysis('Static')

    # Do one analysis for constant axial load
    analyze(1)

    # Define reference moment
    timeSeries('Linear', 2)
    pattern('Plain',2, 2)
    load(2, 0.0, 0.0, 1.0)

    # Compute curvature increment
    dK = maxK / numIncr

    # Use displacement control at node 2 for section analysis
    integrator('DisplacementControl', 2,3,dK,1,dK,dK)

    # Do the section analysis
    analyze(numIncr)
        
    
class LogoController(ViktorController):
    label = 'Logo design'
    parametrization = LogoParametrization

    viktor_convert_entity_field = True
   
    @staticmethod
    def interaction_diagram(params):
        
        #inputs
        analysis_type = params.tab_1.geometry.sectiontype
        concrete_strength = int(params.tab_1.material.concrete)
        steel_strength = int(params.tab_1.material.steel)
        gc = int(params.tab_1.material.gamma_concrete)
        gs = int(params.tab_1.material.gamma_steel)
        long_dia = params.tab_1.longitudinal.diameter
        trans_dia = params.tab_1.transverse.diameter
        b = params.tab_1.geometry.depth
        h = params.tab_1.geometry.width
        cover = params.tab_1.geometry.cover
        d = b-cover
        axial = params.tab_1.force.axialforce
        s = int(params.tab_1.transverse.space)
        nx = params.tab_1.longitudinal.nx
        ny = params.tab_1.longitudinal.ny
        length = params.tab_1.geometry.length
        mb = params.tab_1.shearforce.moment_bottom
        mt = params.tab_1.shearforce.moment_bottom
        ve = params.tab_1.shearforce.earthquake
        v_gqe = params.tab_1.shearforce.combination
        v_gqde = params.tab_1.shearforce.combination_d 
        degree = params.tab_1.geometry.degree
        sectiontype = params.tab_3.axial_frp.sectiontype
        rc = params.tab_3.axial_frp.rc 
        Ef = params.tab_3.axial_frp.youngmodulusofwrap
        ef = params.tab_3.axial_frp.strainofwrap
        eu = params.tab_3.axial_frp.ultimatestrainofwrap
        nf = params.tab_3.axial_frp.numberoffrp
        strip_type = params.tab_3.axial_frp.striptype
        tf = params.tab_3.axial_frp.thicknessofwrap
        sf = params.tab_3.axial_frp.distanceofwraps
        
        if strip_type == "Uncontinuous":
            fi_f = nf*wf*(b+h)*2*tf/(b*h*sf)
        else:
            fi_f = nf*(b+h)*2*tf/(b*h)
            
        if ef <= 0.004:
            ef_1 = ef
        else:
            ef_1 = 0.004
            
        ef_2 = 0.50*eu
        ef = min(ef_1,ef_2)
        
        if sectiontype == "Rectangular Section":
            kappa_a = 1-((b-2*rc)**2 + (h+2*rc)**2)/(3*b*h)
        elif sectiontype == "Circular Section":
            kappa_a = 1
        else:
            kappa_a = b/h
            
  
        f1 = 0.5*fi_f*kappa_a*ef*Ef
        fcc_1 = concrete_strength*(1+2.4*(f1/concrete_strength))
        fcc_2 = 1.2*concrete_strength
        fcc_frp = min(fcc_1,fcc_2) # Compression strength of the concrete to calculate the axial force
        
        concrete_list = [concrete_strength, fcc_frp]
        
        for each in concrete_list:
            #material
            young_modulus_concrete = 5000*sqrt(each)    
            fctk = 0.35*sqrt(each)
            fctd = fctk/gc
            fcd = each/gc
            fyd = steel_strength/gs
            young_modulus_steel = 200000
            s_confinement = max(min(min(b,h)/3,150,long_dia*6),50)
            a = trans_dia*25
            
            diameter_area = long_dia*long_dia*pi/4
            
            x_dir_list = []
            
            if degree == 0:
                d = h - cover
                k = 1
                x_reb = (b-2*cover)/(nx-1)
                increase = (b-2*cover)/(nx-1)
                while k <= nx:
                    if k == 1:
                        x_ = cover
                    elif k == nx:
                        x_ = b - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
            elif degree == 90:
                d = h - cover
                k = 1
                x_reb = (h-2*cover)/(ny-1)
                increase = (h-2*cover)/(ny-1)
                while k <= ny:
                    if k == 1:
                        x_ = cover
                    elif k == ny:
                        x_ = h - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
    
            ecu = 0.003
            esy = fyd/young_modulus_steel
            k1 = 0.82
            index_0 = 15
            index_90 = 15
            k =min(b,h)/index_0
            sensibility_0 = min(b,h)/index_0
            sensibility_90 = min(b,h)/index_90
     
            if degree == 0:
                cb = ecu*(h-cover)/(esy+ecu)
        
                if cb != 0:
                    if cb*k1 > h:
                        k1cb_ = h
                    else:
                        k1cb_ = cb*k1
                total_step = int(ceil((h-cover)/sensibility_0)+2)
                k1cb_list = []
                c_list = []
                deneme_mb_list = []
                c_index = 1
                total_layer_index = 1
                Force_List = []
                Moment_List = []
                while c_index <= total_step:
                    if c_index == 1:
                        c = h*1.3
                    elif c_index == 2:
                        c = h
                    elif c_index == total_step:
                        c = cover
                    else:
                        c = h - sensibility_0   
                        sensibility_0 = sensibility_0+k
                    if c*k1 > h:
                        k1cb = h 
                    else:
                        k1cb = c*k1 
                    k1cb_list.append(k1cb)
            
                    
                    rebar_number_list = []
                    rebar_distance_list = []
                    diameter_list = []
                    area_list = []
                    Le1_list = []
                    ec1_list = []
                    sigma_s_list = []
                    Fs1_list = []
                    Mb_list = []
                    distance_dif_cons = (b-2*cover)/(nx-1)
                    distance_dif = (b-2*cover)/(nx-1)
                    total_layer_index = 1
                    while total_layer_index <= nx:
                        if total_layer_index == 1:
                            rebar_number = ny 
                            rebar_distance = cover
                            area = rebar_number*diameter_area
                           
                        elif total_layer_index > 1 and total_layer_index < nx:
                            rebar_number = 2
                            rebar_distance = cover + distance_dif
                            area = rebar_number*diameter_area
                            distance_dif = distance_dif + distance_dif_cons
                            
                        else:
                            rebar_number = ny
                            rebar_distance = b-cover
                            area = rebar_number*diameter_area
                            
                        Le1 = c - rebar_distance    
                        if (Le1/c)*ecu > 0:
                            if (Le1/c)*ecu > esy:
                                es1 = esy
                            else:
                                es1 = (Le1/c)*ecu
                        elif (Le1/c)*ecu <= 0:
                            if (Le1/c)*ecu < esy*-1:
                                es1 = esy*-1
                            else:
                                es1 = (Le1/c)*ecu
                                
                        if es1*young_modulus_steel > 0:
                            if es1*young_modulus_steel > fyd:
                                sigma = fyd
                            else:
                                sigma = es1*young_modulus_steel
                        elif es1*young_modulus_steel <= 0:
                            if es1*young_modulus_steel <= fyd*-1:
                                sigma = fyd*-1
                            else:
                                sigma = es1*young_modulus_steel
                        Fc = 0.85*fcd*b*k1cb/1000
                              
                        Fs1 = (area*sigma/1000) 
                        Mb_ = Fc*(h/2-k1cb/2)/1000
                        Mb = Fs1*(h/2-rebar_distance)/1000  
                        
                        rebar_number_list.append(rebar_number)
                        rebar_distance_list.append(rebar_distance)
                        diameter_list.append(long_dia)
                        area_list.append(area)
                        Le1_list.append(int(Le1))
                        ec1_list.append(es1)
                        sigma_s_list.append(sigma)
                        Fs1_list.append(Fs1)
                        Mb_list.append(Mb)
                        if c_index == 1:
                            deneme_mb_list.append(Mb)
                
                 
                        #print(Fs1_list)
                        total_layer_index = total_layer_index +1
                    total_Fc = int(ceil(sum(Fs1_list) + Fc))
                    total_Mb = int(ceil(sum(Mb_list) + Mb_))
                    total_area = int(ceil(sum(area_list)))
                    Moment_List.append(total_Mb)
                    Force_List.append(total_Fc)
                    c_list.append(c)
                    c_index = c_index +1 
                    
            elif degree == 90:
                cb = ecu*(b-cover)/(esy+ecu)
                
                if cb != 0:
                    if cb*k1 > b:
                        k1cb_ = b
                    else:
                        k1cb_ = cb*k1
                total_step = int(ceil((b-cover)/sensibility_0)+2)
                k1cb_list = []
                c_list = []
                deneme_mb_list = []
                c_index = 1
                total_layer_index = 1
                Force_List = []
                Moment_List = []
                while c_index <= total_step:
                    if c_index == 1:
                        c = b*1.3
                    elif c_index == 2:
                        c = b
                    elif c_index == total_step:
                        c = cover
                    else:
                        c = b - sensibility_0   
                        sensibility_0 = sensibility_0+k
                    if c*k1 > b:
                        k1cb = b 
                    else:
                        k1cb = c*k1 
                    k1cb_list.append(k1cb)
            
                    
                    rebar_number_list = []
                    rebar_distance_list = []
                    diameter_list = []
                    area_list = []
                    Le1_list = []
                    ec1_list = []
                    sigma_s_list = []
                    Fs1_list = []
                    Mb_list = []
                    distance_dif_cons = (h-2*cover)/(ny-1)
                    distance_dif = (h-2*cover)/(ny-1)
                    total_layer_index = 1
                    while total_layer_index <= ny:
                        if total_layer_index == 1:
                            rebar_number = nx 
                            rebar_distance = cover
                            area = rebar_number*diameter_area
                           
                        elif total_layer_index > 1 and total_layer_index < ny:
                            rebar_number = 2
                            rebar_distance = cover + distance_dif
                            area = rebar_number*diameter_area
                            distance_dif = distance_dif + distance_dif_cons
                            
                        else:
                            rebar_number = nx
                            rebar_distance = h-cover
                            area = rebar_number*diameter_area
                            
                        Le1 = c - rebar_distance    
                        if (Le1/c)*ecu > 0:
                            if (Le1/c)*ecu > esy:
                                es1 = esy
                            else:
                                es1 = (Le1/c)*ecu
                        elif (Le1/c)*ecu <= 0:
                            if (Le1/c)*ecu < esy*-1:
                                es1 = esy*-1
                            else:
                                es1 = (Le1/c)*ecu
                                
                        if es1*young_modulus_steel > 0:
                            if es1*young_modulus_steel > fyd:
                                sigma = fyd
                            else:
                                sigma = es1*young_modulus_steel
                        elif es1*young_modulus_steel <= 0:
                            if es1*young_modulus_steel <= fyd*-1:
                                sigma = fyd*-1
                            else:
                                sigma = es1*young_modulus_steel
                        Fc = 0.85*fcd*h*k1cb/1000
                              
                        Fs1 = (area*sigma/1000) 
                        Mb_ = Fc*(b/2-k1cb/2)/1000
                        Mb = Fs1*(b/2-rebar_distance)/1000  
                        
                        rebar_number_list.append(rebar_number)
                        rebar_distance_list.append(rebar_distance)
                        diameter_list.append(long_dia)
                        area_list.append(area)
                        Le1_list.append(int(Le1))
                        ec1_list.append(es1)
                        sigma_s_list.append(sigma)
                        Fs1_list.append(Fs1)
                        Mb_list.append(Mb)
                        if c_index == 1:
                            deneme_mb_list.append(Mb)
                
                 
                        #print(Fs1_list)
                        total_layer_index = total_layer_index +1
                    total_Fc = int(ceil(sum(Fs1_list) + Fc))
                    total_Mb = int(ceil(sum(Mb_list) + Mb_))
                    total_area = int(ceil(sum(area_list)))
                    Moment_List.append(total_Mb)
                    Force_List.append(total_Fc)
                    c_list.append(c)
                    c_index = c_index +1 
            
            Nc = (b*h*0.85*fcd+fyd*total_area)/1000
            Nt = -1*fyd*total_area/1000
            Moment_List.insert(0, 0)
            Moment_List.insert(total_step+1, 0)
            Force_List.insert(0, Nc)
            Force_List.insert(total_step+1, Nt)
            
            if each == concrete_strength:
                moment_existing = Moment_List
                force_existing = Force_List
            else:
                moment_retrofit = Moment_List
                force_retrofit= Force_List
                
        return moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp
    
    @PlotlyView("Interaction Diagram", duration_guess=1)
    def get_data_view1(self, params, **kwargs):

        analysis_type = params.tab_1.geometry.sectiontype
        moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp = self.interaction_diagram(params)
        
        df_m_n = pd.DataFrame(list(zip(force_existing, moment_existing)), columns =['Axial Force', 'Moment'], dtype = float)
        moment_existing_180 = [x * -1 for x in moment_existing]
        df_m_n_180 = pd.DataFrame(list(zip(force_existing, moment_existing_180)), columns =['Axial Force', 'Moment'], dtype = float)
        moment_retrofit_180 = [x * -1 for x in moment_retrofit]
    
        if analysis_type == "Retrofit":
            fig = go.Figure()
    
            fig.add_trace(go.Scatter(x=moment_existing, y=force_existing,  
                                mode='lines+markers',
                                name='0-Degree_Existing'))
            fig.add_trace(go.Scatter(x=moment_existing_180, y=force_existing,
                                mode='lines+markers',
                                name='180-Degree_Existing'))
            fig.add_trace(go.Scatter(x=moment_retrofit, y=force_retrofit,
                                mode='lines+markers',
                                name='180-Degree_Retrofit'))
            fig.add_trace(go.Scatter(x=moment_retrofit_180, y=force_retrofit,
                                mode='lines+markers',
                                name='180-Degree_Retrofit'))
        else:
            fig = go.Figure()
    
            fig.add_trace(go.Scatter(x=moment_existing, y=force_existing,  
                                mode='lines+markers',
                                name='0-Degree_Existing'))
            fig.add_trace(go.Scatter(x=moment_existing_180, y=force_existing,
                                mode='lines+markers',
                                name='180-Degree_Existing'))
        
        return PlotlyResult(fig.to_json())
    
    @staticmethod
    def calculate_moment_and_curvature(params):
        
        moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp = LogoController.interaction_diagram(params)
        concrete_strength = int(params.tab_1.material.concrete)
        steel_strength = int(params.tab_1.material.steel)
        gc = int(params.tab_1.material.gamma_concrete)
        gs = int(params.tab_1.material.gamma_steel)
        long_dia = params.tab_1.longitudinal.diameter
        trans_dia = params.tab_1.transverse.diameter
        b = params.tab_1.geometry.depth
        h = params.tab_1.geometry.width
        depth = params.tab_1.geometry.depth
        width = params.tab_1.geometry.width
        cover = params.tab_1.geometry.cover
        d = b-cover
        axialLoad = params.tab_1.force.axialforce
        s = int(params.tab_1.transverse.space)
        nx = params.tab_1.longitudinal.nx
        ny = params.tab_1.longitudinal.ny
        length = params.tab_1.geometry.length
        mb = params.tab_1.shearforce.moment_bottom
        mt = params.tab_1.shearforce.moment_bottom
        ve = params.tab_1.shearforce.earthquake
        v_gqe = params.tab_1.shearforce.combination
        v_gqde = params.tab_1.shearforce.combination_d 
        degree = params.tab_1.geometry.degree
        n_leg_x = params.tab_1.transverse.nlegx
        n_leg_y = params.tab_1.transverse.nlegy
        ax1 = params.tab_2.concreteproperties.ax
        ay1 = params.tab_2.concreteproperties.ay
        coefficient = str(params.tab_1.material.coefficient)
        s = params.tab_1.transverse.space
        eu = params.tab_3.axial_frp.ultimatestrainofwrap
        analysis_type = params.tab_1.geometry.sectiontype
        confinement_type = params.tab_1.transverse.confinementtype
        
        #material
        young_modulus_concrete = 5000*sqrt(concrete_strength)    
        fctk = 0.35*sqrt(concrete_strength)
        fctd = fctk/gc
        fcd = concrete_strength/gc
        fyd = steel_strength/gs
        young_modulus_steel = 200000
        s_confinement = max(min(min(b,h)/3,150,long_dia*6),50)
        a = trans_dia*25
        
        ef_1 = 0.50*eu            
        ef_2 = 0.01
        ef = min(ef_1,ef_2)
        
        ecc_frp = 0.002*(1+15*((f1/fcd)**0.75))
        fcc_retrofit = fcc_frp
        
        diameter_area = long_dia*long_dia*pi/4

        design_list = ["Existing", "Retrofit"]
        
        for each in design_list:
            design_type = each
            fcc, eco, esp, ecu_max, fcu_max, df_conf, eco_max, ecu_maximum, ecc = concrete_func(design_type, ecc_frp, fcc_retrofit, ax1, ay1, depth,width, cover, long_dia, trans_dia, nx, ny, n_leg_x, n_leg_y, coefficient, concrete_strength, steel_strength, s)
            x_dir_list = []
    
            if degree == 0:
                d = h - cover
                k = 1
                x_reb = (b-2*cover)/(nx-1)
                increase = (b-2*cover)/(nx-1)
                while k <= nx:
                    if k == 1:
                        x_ = cover
                    elif k == nx:
                        x_ = b - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
            elif degree == 90:
                d = h - cover
                k = 1
                x_reb = (h-2*cover)/(ny-1)
                increase = (h-2*cover)/(ny-1)
                while k <= ny:
                    if k == 1:
                        x_ = cover
                    elif k == ny:
                        x_ = h - cover
                    else:
                        x_ = cover + increase
                        increase = increase + x_reb
                    x_dir_list.append(int(x_))
                    k = k + 1
            
            ecu = 0.002
            esy = fyd/young_modulus_steel
            
            if degree == 0:
                y1 = int(b/2.0)
                z1 = int(h/2.0)
            
                length_le1 = len(x_dir_list)          
                
            elif degree == 90:
                y1 = int(h/2.0)
                z1 = int(b/2.0)
                length_le1 = len(x_dir_list)          
    
            wipe()
            print("Start MomentCurvature.py example")
            
            concrete_output = "concrete_strain.txt"
            steel_output = "steel_strain.txt"
            moment_curvature_output = "moment_curvature.txt"   
            name = "name"
            
            # ------------------------------ Moment Curvature ------------------------------
            # Define model builder
            # --------------------
            model('basic','-ndm',2,'-ndf',3)                   
            secTag = 1        
            
            # Cover concrete (unconfined)
            # uniaxialMaterial('Concrete04',2, -fcd,  -eco,  -esp,  young_modulus_concrete, 0.0, 0.0, 0,1)
            # uniaxialMaterial('Concrete01',1, -fcc, -eco_max,  -fcu_max,  -ecu_maximum)
            uniaxialMaterial('Concrete02',1, -fcu_max, -ecu,  -0.2*fcu_max,  -0.02, 0.1, 0.1*fcu_max, (0.1*fcu_max)/ecu)
            # uniaxialMaterial('Concrete01',1, -fcc, -eco_max,  -fcu_max,  -0.02)
            # uniaxialMaterial('Concrete04',1, float(fcc),  ecc,  -0.02,  Ec)
            
            # # Cover concrete (unconfined)
            # uniaxialMaterial('Concrete04',2, float(fc0),  -0.002,  -0.004,  Ec)
            uniaxialMaterial('Concrete02',2, -fcd, -ecu,  -0.2*fcd,  -0.02, 0.1, 0.1*fcd, (0.1*fcd)/ecu)
            # uniaxialMaterial('Concrete01',2, -fcd,  -ecu,  0.0,  -esp)
            
            # STEEL
            # Reinforcing steel 
            by = 0.01
            
            #                        tag  fy E0    b
            uniaxialMaterial('Steel01', 3, int(fyd), young_modulus_steel, by)
            
            # Define cross-section for nonlinear columns
            # ------------------------------------------
            
            if degree == 0:
                width = h
                depth = b
                cover = cover
                number_of_layer = length_le1
                n_top = ny
                n_bot = ny
                n_int = 2
                
                
                b1 = width/2 - cover
                b2 = (width/2 - cover)*-1
                h1 = depth/2 - cover
                h2 = (depth/2 - cover)*-1
                
                # some variables derived from the parameters
                y1 = depth/2.0
                z1 = width/2.0
                total_y = depth - 2*cover
                total_y_layer = total_y/(number_of_layer-1)
                total_y_layer_step = total_y/(number_of_layer-1)
                
                section('Fiber', 1)
                
                # Create the concrete core fibers
                if confinement_type == "Confined" or confinement_type == "Unconfined" and design_type == "Retrofit":
                    patch('rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                    
                    
                    # Create the concrete cover fibers (top, bottom, left, right)
                    patch('rect',2,20,1 ,-y1, z1-cover, y1, z1)
                    patch('rect',2,20,1 ,-y1, -z1, y1, cover-z1)
                    patch('rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover)
                    patch('rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover)

                
                    top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                    bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                    
                    fib_sec_1 = [['section', 'Fiber', 1],
                    ['patch', 'rect',2,20,1 ,-y1, z1-cover, y1, z1],
                    ['patch', 'rect',2,20,1 ,-y1, -z1, y1, cover-z1],
                    ['patch', 'rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover],
                    ['patch', 'rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover],
                    ['patch', 'rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                    top,
                    bottom]
                    
                elif confinement_type == "Unconfined" and design_type == "Existing":
                    patch('rect',2,10,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                
                    top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                    bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                    
                    fib_sec_1 = [['section', 'Fiber', 1],
                    ['patch', 'rect',2,10,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                    top,
                    bottom]
                
                elif confinement_type == "Confined" and design_type == "Existing":
                    patch('rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                    
                    
                    # Create the concrete cover fibers (top, bottom, left, right)
                    patch('rect',2,20,1 ,-y1, z1-cover, y1, z1)
                    patch('rect',2,20,1 ,-y1, -z1, y1, cover-z1)
                    patch('rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover)
                    patch('rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover)

                
                    top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                    bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                    
                    fib_sec_1 = [['section', 'Fiber', 1],
                    ['patch', 'rect',2,20,1 ,-y1, z1-cover, y1, z1],
                    ['patch', 'rect',2,20,1 ,-y1, -z1, y1, cover-z1],
                    ['patch', 'rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover],
                    ['patch', 'rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover],
                    ['patch', 'rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                    top,
                    bottom]
                    
                # plt.savefig('fibsec_rc.png')
                
                # # Create the reinforcing fibers (left, middle, right)
                
                layer('straight', 3, n_top, diameter_area, y1-cover, cover-z1, y1-cover, z1-cover)
                layer('straight', 3, n_bot, diameter_area, cover-y1, cover-z1, cover-y1, z1-cover)
                
                total_int_layer = number_of_layer-2
                int_layer = 1
                while int_layer <= total_int_layer:
                
                    layer('straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia))
                    int_layer_def = ['layer','straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia)]
                    fib_sec_1.append(int_layer_def)
                    total_y_layer = total_y_layer + total_y_layer_step
                    int_layer = int_layer +1
                    
                # d -- from cover to rebar
                d = depth-cover
            
            elif degree == 90:
                width = b
                depth = h
                cover = cover
                number_of_layer = length_le1
                n_top = nx
                n_bot = nx
                n_int = 2
                
                
                b1 = width/2 - cover
                b2 = (width/2 - cover)*-1
                h1 = depth/2 - cover
                h2 = (depth/2 - cover)*-1
                
                # some variables derived from the parameters
                y1 = depth/2.0
                z1 = width/2.0
                total_y = depth - 2*cover
                total_y_layer = total_y/(number_of_layer-1)
                total_y_layer_step = total_y/(number_of_layer-1)
                
                section('Fiber', 1)
                
                # Create the concrete core fibers
                patch('rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover)
                
                
                # Create the concrete cover fibers (top, bottom, left, right)
                patch('rect',2,20,1 ,-y1, z1-cover, y1, z1)
                patch('rect',2,20,1 ,-y1, -z1, y1, cover-z1)
                patch('rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover)
                patch('rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover)
                
                top = ['layer','straight', 3, n_top, diameter_area, y1-cover-long_dia, cover-z1+long_dia, y1-cover-long_dia, z1-cover-long_dia]
                bottom = ['layer','straight', 3, n_bot, diameter_area, cover-y1+long_dia, cover-z1+long_dia, cover-y1+long_dia, z1-cover-long_dia]
                
                fib_sec_1 = [['section', 'Fiber', 1],
                ['patch', 'rect',2,20,1 ,-y1, z1-cover, y1, z1],
                ['patch', 'rect',2,20,1 ,-y1, -z1, y1, cover-z1],
                ['patch', 'rect',2,2,1 ,-y1, cover-z1, cover-y1, z1-cover],
                ['patch', 'rect',2,2,1 , y1-cover, cover-z1, y1, z1-cover],
                ['patch', 'rect',1,20,1 ,cover-y1, cover-z1, y1-cover, z1-cover],
                top,
                bottom]
                    
                
                # plt.savefig('fibsec_rc.png')
                
                # # Create the reinforcing fibers (left, middle, right)
                
                layer('straight', 3, n_top, diameter_area, y1-cover, cover-z1, y1-cover, z1-cover)
                layer('straight', 3, n_bot, diameter_area, cover-y1, cover-z1, cover-y1, z1-cover)
                
                total_int_layer = number_of_layer-2
                int_layer = 1
                while int_layer <= total_int_layer:
                
                    layer('straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia))
                    int_layer_def = ['layer','straight', 3, int(n_int), int(diameter_area), int(y1-cover-total_y_layer), int(cover-z1+long_dia), int(y1-cover-total_y_layer), int(z1-cover-long_dia)]
                    fib_sec_1.append(int_layer_def)
                    total_y_layer = total_y_layer + total_y_layer_step
                    int_layer = int_layer +1
                    
                # d -- from cover to rebar
                d = width-cover
                
            matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
            m_c = opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor) #Figure of the section
            
            plt.axis('equal')    
            
            # Estimate yield curvature
            # (Assuming no axial load and only top and bottom steel)
            # d -- from cover to rebar
            
            # steel yield strain
            epsy = fyd/young_modulus_steel
            Ky = epsy/(0.7*d)
            
            # Print estimate to standard output
            print("Estimated yield curvature: ", Ky)
               
            # Target ductility for analysis
            mu = 30
            
            # Number of analysis increments
            numIncr = 500
            
            # Call the section analysis procedure
            
            results = open('results.out','a+')
            
            # u = nodeDisp(2,3)
            # if abs(u-0.00190476190476190541)<1e-12:
            #     results.write('PASSED : MomentCurvature.py\n');
            #     print("Passed!")
            # else:
            #     results.write('FAILED : MomentCurvature.py\n');
            #     print("Failed!")
            
            results.close()
            
            print("==========================")
            
            axialDesign = -1000*axialLoad
            
            MomentCurvature(name, concrete_output,steel_output,moment_curvature_output, secTag, b1, b2, h1, h2, axialDesign, Ky*mu, numIncr)
            
            
            # Reading of Moment Curvature Results
            with open(moment_curvature_output) as f:
                coords = f.read().split()
            
            # Splitting data as Moment & Curvature
            moment = coords[0::2]
            curvature = coords[1::2]
            
            moment = moment[:-1]
            curvature = curvature[:-1]
            
            df = pd.DataFrame(list(zip(curvature, moment)), columns =['Curvature', 'Moment'], dtype = float)
            df['Curvature'] = 1000*df['Curvature']
            df['Moment'] = df['Moment']/1000000
            df.plot(kind='line',x='Curvature',y='Moment',color='red')
            moment_list = df['Moment'].tolist()
            curvature_list = df["Curvature"].tolist()
            
            if each == "Existing":
                moment_existing = moment_list
                curvature_existing = curvature_list
                df_existing = pd.DataFrame(list(zip(curvature_existing, moment_existing)), columns =['Curvature', 'Moment'], dtype = float)
            elif each == "Retrofit":
                moment_retrofit = moment_list
                curvature_retrofit= curvature_list
                df_retrofit = pd.DataFrame(list(zip(curvature_retrofit, moment_retrofit)), columns =['Curvature', 'Moment'], dtype = float)
            
            
        
        def convert_design_m_c(df):
            return df.to_csv().encode('utf-8')

        csv_design = convert_design_m_c(df)
        document_code_mc = "Moment - Curvature (0 Degree) - Press to Download"
        csv_file_mc = "m_n_interaction_0.csv"
        
            
        return moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing
    
    @PlotlyView("Moment Curvature", duration_guess=1)
    def get_plotly_view(self, params, **kwargs):
        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing = self.calculate_moment_and_curvature(params)
        
        fig = go.Figure()
        # fig = px.line(x=Moment_List, y=Force_List, labels={'x':'Moment', 'y':'Axial Force'})
        # Add traces
        # fig.add_trace(go.Scatter(x=moment_list, y=curvature_list,  
        #                     mode='lines+markers',
        #                     name='0-Degree'))
        # # fig.add_trace(go.Scatter(x=curvature_list, y=moment_list,
        # #                     mode='lines+markers',
        # #                     name='180-Degree'))
        

        fig = go.Figure()

        fig.add_trace(go.Scatter(x=curvature_existing, y=moment_existing,  
                            mode='lines+markers',
                            name='0-Degree_Existing'))
        fig.add_trace(go.Scatter(x=curvature_retrofit, y=moment_retrofit,
                            mode='lines+markers',
                            name='0-Degree_Retrofit'))
        
        #fig = px.line(x=curvature_list, y=moment_list, labels={'x':'Curvature', 'y':'Moment'})
        
        return PlotlyResult(fig.to_json())
    

    @GeometryView('3D', duration_guess=1)
    def get_3d_view(self, params, **kwargs):

        depth = params.tab_1.geometry.depth
        width = params.tab_1.geometry.width
        height = params.tab_1.geometry.length

        beam = SquareBeam(depth, width, height)

        # axis_system_location = Point(2, 2, 2)
        # axis_system = CartesianAxes(axis_system_location, axis_length=0.5, axis_diameter=0.5)

        geometries = [beam]

        return GeometryResult(geometries)
    
    def download_mass_breakdown_existing(self, params, **kwargs):

        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing = self.calculate_moment_and_curvature(params)

        # df = pd.DataFrame(moment_list)
        
        return DownloadResult(df_existing.to_csv(), 'moment_curvature_existing.csv')
    
    def download_mass_breakdown_retrofit(self, params, **kwargs):

        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum, df_retrofit, df_existing = self.calculate_moment_and_curvature(params)

        # df = pd.DataFrame(moment_list)
        
        return DownloadResult(df_retrofit.to_csv(), 'moment_curvature_retrofit.csv')
    
    @DataView('Shear Capacity', duration_guess=1)
    def get_data_view(self, params, **kwargs):

        #inputs
        moment_existing, force_existing, moment_retrofit, force_retrofit, f1, fcc_frp = self.interaction_diagram(params)
        moment_existing, curvature_existing, moment_retrofit, curvature_retrofit, df, ecc, ecu_maximum = self.calculate_moment_and_curvature(params)
        concrete_strength = int(params.tab_1.material.concrete)
        steel_strength = int(params.tab_1.material.steel)
        gc = int(params.tab_1.material.gamma_concrete)
        gs = int(params.tab_1.material.gamma_steel)
        long_dia = params.tab_1.longitudinal.diameter
        trans_dia = params.tab_1.transverse.diameter
        depth = params.tab_1.geometry.depth
        width = params.tab_1.geometry.width
        cover = params.tab_1.geometry.cover
        d = depth-cover
        axial = params.tab_1.force.axialforce
        s = int(params.tab_1.transverse.space)
        n_leg_x = params.tab_1.transverse.nlegx
        n_leg_y = params.tab_1.transverse.nlegy
        length = params.tab_1.geometry.length
        mb = params.tab_1.shearforce.moment_bottom
        mt = params.tab_1.shearforce.moment_bottom
        ve = params.tab_1.shearforce.earthquake
        v_gqe = params.tab_1.shearforce.combination
        v_gqde = params.tab_1.shearforce.combination_d 
        degree = params.tab_1.geometry.degree
        nf = params.tab_3.frp.numberoffrp
        tf = params.tab_3.frp.thicknessofwrap
        wf = params.tab_3.frp.widthofwrap
        Ef = params.tab_3.frp.youngmodulusofwrap
        ef = params.tab_3.frp.strainofwrap
        d_wrap = params.tab_3.frp.effectivedepth
        sf = params.tab_3.frp.distanceofwraps
        strip_type = params.tab_3.frp.striptype
        eu = params.tab_3.frp.ultimatestrainofwrap
    
        #material
        young_modulus_concrete = 5000*sqrt(concrete_strength)    
        fctk = 0.35*sqrt(concrete_strength)
        fctd = fctk/gc
        fcd = concrete_strength/gc
        fyd = steel_strength/gs
        young_modulus_steel = 200000
        s_confinement = max(min(min(depth,width)/3,150,long_dia*6),50)
        a = trans_dia*25
        
        if axial < 0:
            gamma = -0.3
        else:
            gamma = 0.07

        #geometry
        
        axial_limit = depth*width*0.2*concrete_strength/1000
        axial_limit_ts500 = depth*width*0.05*concrete_strength/1000
        depth_clear = depth-2*cover-trans_dia
        width_clear = width-2*cover-trans_dia
        area_clear = depth_clear*width_clear
        area = depth*width
        trans_area_x = pi*trans_dia*trans_dia*n_leg_x/4 
        trans_area_y = pi*trans_dia*trans_dia*n_leg_y/4 
        s_mid = min(min(depth,width)/2,200)
        
        #shear force
        ve_capacity = (mb+mt)/length
        
        ve_final = min(ve_capacity, v_gqde)
        
        if axial < axial_limit_ts500 and ve > v_gqe/2:
            vc = 0
        else:
            vc = ((0.65*fctd*d*depth)*(1+gamma*axial*1000/(depth*d))*0.001)*0.8
        vw = trans_area_x*fyd*d/(s*1000)
        vr = vc + vw
        
        v_tbdy2018 = 0.85*depth*d*sqrt(fcd)*0.001
        
        if vr > ve_final and v_tbdy2018 > ve_final:
            shear_check = "TS500 and TBDY 2018 check are ok!"
        elif vr > ve_final and v_tbdy2018 < ve_final:
            shear_check = "TS500 check is ok - TBDY 2018 check is not ok!"
        elif vr < ve_final and v_tbdy2018 > ve_final:
            shear_check = "TS500 check is not ok - TBDY 2018 check is ok!"   
        elif vr < ve_final and v_tbdy2018 < ve_final:
            shear_check = "TS500 and TBDY 2018 check are not ok!" 

        
        if axial > axial_limit:
            con_reb_area_x = max(0.30*s_confinement*depth_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*depth_clear*(fcd/fyd))
            con_reb_area_y = max(0.30*s_confinement*width_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*width_clear*(fcd/fyd))
        else:
            con_reb_area_x = max(0.30*s_confinement*depth_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*depth_clear*(fcd/fyd))*2/3
            con_reb_area_y = max(0.30*s_confinement*width_clear*((area/area_clear)-1)*(fcd/fyd),0.075*s*width_clear*(fcd/fyd))*2/3
           
        if trans_area_x < con_reb_area_x:
            warning = "Transverse spacing is not enough for confinement zone!"
        else:
            warning = "Transverse spacing is enough for confinement zone."
        
        if ef <= 0.004:
            ef_1 = ef
        else:
            ef_1 = 0.004
            
        ef_2 = 0.50*eu
        ef = min(ef_1,ef_2)
         
        if strip_type == "Uncontinuous":
            sf = min(wf+d/4, sf)
            v_wrap = 2*nf*tf*wf*Ef*ef*d_wrap/sf/1000
        else:
            v_wrap = 2*nf*tf*Ef*ef*d_wrap/1000
        
        data = DataGroup(
            spacing_conf = DataItem(warning, s, suffix='mm', number_of_decimals=2),
            spacing_mid = DataItem("Transverse space is: ", s_mid, suffix='mm', number_of_decimals=2),
            tbdy_2018 = DataItem("TBDY2018 - Capacity:", v_tbdy2018, suffix='kN',number_of_decimals=2),
            # ef_final = DataItem("Unit Strain of FRP:", ef, number_of_decimals=3),
            # ecc_final = DataItem("Unit Strain of FRP:", ecc, number_of_decimals=3),
            # ecu_maximum1 = DataItem("Ecu Max:", ecu_maximum, number_of_decimals=3),
            # f1 = DataItem("f1:", f1, number_of_decimals=3),
            # fcc_frp = DataItem("fcc_frp:", fcc_frp, number_of_decimals=3),
            shearcapacity_stirrups = DataItem("Shear Force Capacity (Vw): ", vw, suffix='kN', number_of_decimals=2),
            shearcapacity_concrete = DataItem("Shear Force Capacity (Vc): ", vc, suffix='kN', number_of_decimals=2),
            frp_capacity = DataItem("Shear Capacity (Vf): ", v_wrap, suffix='kN', number_of_decimals=2),
            shearforce_capacity = DataItem("Total Shear Force Capacity without FRP (Vr): ", vr, suffix='kN', number_of_decimals=2),
            shearforce_capacity_frp = DataItem("Total Shear Force Capacity with FRP (Vf): ", v_wrap+vr, suffix='kN', number_of_decimals=2),
            
            # capacity_check = DataItem("Shear Capacity Check:", shear_check)

        )
        
        return DataResult(data)
    
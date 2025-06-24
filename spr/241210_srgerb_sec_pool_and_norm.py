from opentrons import protocol_api
import numpy as np

#Script originally written by Basile Wicky <bwicky@uw.edu>
#Maintained by Stacey Gerben <srgerb@uw.edu> 
#Let me know if there are bugs

metadata = {
    'protocolName': '241210_srgerb_sec_pool_and_norm',
    'author': 'srgerb',
    'apiLevel':'2.5',
    'description': 'Pool SEC fractions and normalize concentrations (optional)',
}

#======================
# TRANSFER DEFINITIONS
#======================

instrument = 'hplc' # {akta or hplc or cwby} instrument used for SEC -- determines the plate-type of the fractions. Use cwby if normalizing pre-pooled fractions.
normalize = True # whether or not to normalize concentrations of pooled fractions.
mixing = True # set to False to overwrite the automatic assignemnt of this variable.
debug = False # use for simulation purposes only.

"""
Format of transfer definitions (list of lists):
[
[['frac_1', 'frac_2', ...], [frac_1_vol, frac_2_vol, ...], 'destination_plate', 'destination_well', buffer_volume],
]

Format of frac_x strings:
'akta':'1.A.1'
'hplc':'1-P1-A1'
'cwby':'1_A1'
"""

transfers = np.array([
[ ['1-P1-A2',], [87.64405249071999], 1, 'A1', 1412.356],
[ ['1-P1-A8',], [159.64000000000013], 1, 'A3', 663.099],
[ ['1-P1-A12',], [85.54298954522594], 1, 'A4', 1414.457],
[ ['1-P1-A15',], [116.54264898022414], 1, 'A5', 1383.457],
[ ['1-P1-A18',], [159.64000000000013], 1, 'A6', 321.597],
[ ['1-P1-B22',], [35.827016957721305], 1, 'A7', 1464.173],
[ ['1-P1-B14',], [111.47596846485877], 1, 'A8', 1388.524],
[ ['1-P1-B8',], [137.31639376457392], 1, 'A9', 1362.684],
[ ['1-P1-C3',], [71.37355896273047], 1, 'A10', 1428.626],
[ ['1-P1-C9',], [160.15999999999997], 1, 'A11', 104.464],
[ ['1-P1-C15',], [47.987700630645335], 1, 'A12', 1452.012],
[ ['1-P1-C21',], [111.3025951180848], 1, 'B1', 1388.697],
[ ['1-P1-D20',], [60.068176825793586], 1, 'B2', 1439.932],
[ ['1-P1-D10',], [37.36994198084668], 1, 'B3', 1462.630],
[ ['1-P1-D4',], [196.9776995099922], 1, 'B4', 1303.022],
[ ['1-P1-E3',], [84.71709347408819], 1, 'B5', 1415.283],
[ ['1-P1-E7',], [199.98114794173472], 1, 'B6', 1300.019],
[ ['1-P1-E13',], [160.15999999999997], 1, 'B7', 475.762],
[ ['1-P1-E19',], [160.15999999999997], 1, 'B8', 830.359],
[ ['1-P1-F21',], [48.45113383918978], 1, 'B9', 1451.549],
[ ['1-P1-F11',], [83.02206947276426], 1, 'B10', 1416.978],
[ ['1-P1-F3',], [178.55174488327202], 1, 'B11', 1321.448],
[ ['1-P1-G8',], [56.38745492300265], 1, 'B12', 1443.613],
[ ['1-P1-G17',], [92.83570663211546], 1, 'C1', 1407.164],
[ ['1-P1-H24',], [116.08047904625069], 1, 'C2', 1383.920],
[ ['1-P1-H16',], [160.15999999999997], 1, 'C3', 1009.991],
[ ['1-P1-H8',], [33.96268521955988], 1, 'C4', 1466.037],
[ ['1-P1-H2',], [101.51211670723067], 1, 'C5', 1398.488],
[ ['1-P1-I6',], [73.27048226760027], 1, 'C6', 1426.730],
[ ['1-P1-I16',], [94.60661379346574], 1, 'C7', 1405.393],
[ ['1-P1-J24',], [129.1767373401992], 1, 'C8', 1370.823],
[ ['1-P1-J17',], [138.11427233996574], 1, 'C9', 1361.886],
[ ['1-P1-J10',], [101.35906023451336], 1, 'C10', 1398.641],
[ ['1-P1-J1',], [160.16000000000014], 1, 'C11', 429.316],
[ ['1-P1-K7',], [34.457635376006586], 1, 'C12', 1465.542],
[ ['1-P1-K13',], [160.15999999999997], 1, 'D1', 528.892],
[ ['1-P1-K22',], [88.38030781409975], 1, 'D2', 1411.620],
[ ['1-P1-L22',], [159.920328601351], 1, 'D3', 1340.080],
[ ['1-P1-L16',], [63.27415624889028], 1, 'D4', 1436.726],
[ ['1-P1-L10',], [61.44593502942983], 1, 'D5', 1438.554],
[ ['1-P1-M2',], [110.12378598837157], 1, 'D6', 1389.876],
[ ['1-P1-M7',], [160.16000000000014], 1, 'D7', 73.235],
[ ['1-P1-M13',], [45.67349425406983], 1, 'D8', 1454.327],
[ ['1-P1-M22',], [32.67795634505984], 1, 'D9', 1467.322],
[ ['1-P1-N20',], [160.16000000000014], 1, 'D10', 1030.753],
[ ['1-P1-N15',], [162.71026477990404], 1, 'D11', 1337.290],
[ ['1-P1-N10',], [160.16000000000014], 1, 'D12', 6.706],

], dtype=object)

def run(protocol: protocol_api.ProtocolContext):

    protocol.set_rail_lights(True)

    # Global defintions.
    z_off = 0.5 # mm -- tip offset from bottom of well.
    z_in = 3 # mm -- tip insertion below the liquid surface.

    # Labware definitions.
    s300 = protocol.load_instrument('p300_single_gen2', 'left')
    s300.flow_rate.aspirate = 300 # uL/s
    s300.flow_rate.dispense = 300 # uL/s
    s300.flow_rate.blow_out = 300 # uL/s
    
    if debug: # use debug plate definitions for simulation purposes only.
        plate_types = {
            'akta':'usascientific_96_wellplate_2.4ml_deep',
            'hplc':'corning_384_wellplate_112ul_flat',
            'cwby':'corning_96_wellplate_360ul_flat',
        }  
    
    else:
        plate_types = {
            'akta':'abgene_96_wellplate_2200ul',
            'hplc':'greiner_384_wellplate_250ul',
            'cwby':'nunc_96_wellplate_450ul',
        } 
    
    # Determine the number of destination plates.
    des_plates = np.unique(transfers[:, 2])
  
    # Determine the maximum volume of any well in the pooled plate
    if normalize:
        max_pool_vol = np.max([np.sum(t[1]) + t[4] for t in transfers])

    else: 
        max_pool_vol = np.max([np.sum(t[1]) for t in transfers])

    if max_pool_vol > 2000:
        protocol.pause('CAUTION!!! Some pooled wells will overflow!!!')

    # Determine the number of source plates.
    if instrument == 'akta':
        src_plates = list(sorted(np.unique([x.split('.')[0] for x in np.hstack(transfers[:, 0])])))
    
    elif instrument == 'hplc':
        src_plates = list(sorted(np.unique(['-'.join(x.split('-')[:2]) for x in np.hstack(transfers[:, 0])])))

    elif instrument == 'cwby':
        src_plates = list(sorted(np.unique([x.split('_')[0] for x in np.hstack(transfers[:, 0])])))

    else:
        print('Instrument not recognized.')
      
    # Assign plate-type and plate locations to deck.
    plate2deck = {0:10, 1:11, 2:7, 3:8, 4:4, 5:5}
    frac_ptype = plate_types[instrument]
    frac_plates = {}
    for i, p_name in enumerate(src_plates):  
        frac_plates[p_name] = protocol.load_labware(frac_ptype, plate2deck[i], label=f'Fractions {p_name}\n{frac_ptype}')
       
    tips = {}
    pool_ptype = plate_types['akta'] if max_pool_vol > 405 else plate_types['cwby']
    pool_plates = {}
    for i, p_name in enumerate(des_plates):
        tips[p_name] = protocol.load_labware('opentrons_96_tiprack_300ul', 9-(i*3), label=f'Tips {p_name}')
        pool_plates[p_name] = protocol.load_labware(pool_ptype, i+1, label=f'Destination {p_name}\n{pool_ptype}')
      

    def V2H(vol, z_insert, z_offset, plate):
        '''
        For 96-deepwell plate (AB0932); 39 x 8.4 x 8.4 mm
        For 384-deepwell plate (Greiner 781270); 19.3 x 3.8 x 3.8 mm
        '''
        h_factor = {
                'usascientific_96_wellplate_2.4ml_deep' :0.0177273,
                'corning_384_wellplate_112ul_flat'      :0.0804167,
                'corning_96_wellplate_360ul_flat'       :0.0,
                'abgene_96_wellplate_2200ul'            :0.0177273,
                'greiner_384_wellplate_250ul'           :0.0804167,
                'nunc_96_wellplate_450ul'               :0.0,

        }

        h = h_factor[plate] * vol # height of the liquid.
        h_corr = h - z_insert # height corrected for tip insertion.
        h_actual = h_corr if h_corr > z_offset else z_offset # ensure that min height is z_offset.
        return h_actual

    # Add buffer if normalizing concentrations.
    if normalize:
        buffer = protocol.load_labware('agilent_1_reservoir_290ml', 3, label='Buffer')
        s300.pick_up_tip(tips[des_plates[0]]['A1'])
        s300.aspirate(300, buffer['A1'].bottom(4))

        for t in transfers:
            _, _, des_p, des_w, buffer_vol = t

            if buffer_vol == 0:
                pass

            else:
            
                if buffer_vol <= s300.current_volume:
                    s300.dispense(buffer_vol, pool_plates[des_p][des_w].bottom(z_off))

                else:

                    if buffer_vol <= 300: 
                        s300.aspirate(300-s300.current_volume, buffer['A1'].bottom(4))
                        s300.dispense(buffer_vol, pool_plates[des_p][des_w].bottom(z_off))

                    else:
                        full_transfers = int(buffer_vol // 300)
                        remaining_vol = buffer_vol - (full_transfers * 300)
                        added_vol = 0
                        for _ in range(full_transfers):
                            s300.aspirate(300-s300.current_volume, buffer['A1'].bottom(4))
                            added_vol += 300
                            s300.dispense(300, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
                        
                        s300.aspirate(300, buffer['A1'].bottom(4))
                        added_vol += remaining_vol
                        s300.dispense(remaining_vol, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
        
        s300.blow_out(buffer['A1'])
        s300.return_tip()
    
    # Pool.   
    for t in transfers:
        fractions, frac_vol, des_p, des_w, buffer_vol = t
        tot_vol = np.sum(frac_vol) + buffer_vol
        s300.pick_up_tip(tips[des_p][des_w])
        added_vol = buffer_vol

        for frac, vol in list(zip(fractions, frac_vol)):

            if instrument == 'hplc':
                src_p, src_w = '-'.join(frac.split('-')[:2]), frac.split('-')[2]

            elif instrument == 'akta':
                src_p, src_w = frac.split('.')[0], ''.join(frac.split('.')[1:])

            elif instrument == 'cwby':
                src_p, src_w = frac.split('_')
            
            else:
                print('Instrument not recognized.')

            
            if frac_ptype == 'greiner_384_wellplate_250ul':
                # Step draw to avoid overflowing the well.
                # Only setup for taking everything from the well.
                #changed 10-29-24 to take less than everything from the well
                if vol < 50 and buffer_vol > 0:
                    s300.aspirate(vol, frac_plates[src_p][src_w].bottom(8.3))
                elif vol < 50:
                    s300.aspirate(vol, frac_plates[src_p][src_w].bottom(z_off))
                elif vol > 50:
                    s300.aspirate(50, frac_plates[src_p][src_w].bottom(8.3)) 
                    s300.aspirate(vol-50, frac_plates[src_p][src_w].bottom(z_off))
                added_vol += vol
                s300.dispense(250, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
                
            else:
                full_transfers = int(vol // 300)
                remaining_vol = vol - (full_transfers * 300)
                for _ in range(full_transfers):
                    s300.aspirate(300, frac_plates[src_p][src_w].bottom(z_off))
                    added_vol += 300
                    s300.dispense(300, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
           
                s300.aspirate(remaining_vol, frac_plates[src_p][src_w].bottom(z_off))
                added_vol += remaining_vol
                s300.dispense(remaining_vol, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))


                                                           
        if ((len(fractions) > 1) | ((buffer_vol != 0) & normalize)) & mixing: # mix if multiple fractions pooled or buffer was added.
            s300.mix(
                (lambda x: 3 if x <= 600 else (1.5 * tot_vol) // 300)(tot_vol), 
                (lambda x: x / 2 if x <= 600 else 300)(tot_vol), 
                pool_plates[des_p][des_w].bottom(V2H(tot_vol, z_in, z_off, pool_ptype))
            ) 

        s300.drop_tip()
    
    protocol.set_rail_lights(False)
from opentrons import protocol_api

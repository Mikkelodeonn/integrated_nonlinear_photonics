import traceback
import emodeconnection as emode
 
def simulate_waveguide(em, 
                       wavelength = 1550, 
                       core_height = 136, 
                       core_width = 1200, 
                       clad_width = 1200, 
                       clad_height = 1200, 
                       dx = 10, 
                       dy = 2, 
                       nr_modes = 1, 
                       max_index = 0, 
                       BC = 'TE', 
                       field_string = 'Ex', 
                       plot_bool = False, 
                       return_fields=False):
    """
    Simulates an InGaP waveguide using the emodeconnection library.
 
    Parameters:
    - em: An instance of the emodeconnection class.
    - wavelength: Wavelength of light in nm.
    - core_height: Height of the waveguide core in nm.
    - core_width: Width of the waveguide core in nm.
    - clad_width: Width of the cladding in nm.
    - clad_height: Height of the cladding in nm.
    - dx: Resolution in x-direction in nm.
    - dy: Resolution in y-direction in nm.
    - nr_modes: Number of modes to compute.
    - max_index: Maximum index for mode computation.
    - BC: Boundary condition ('TE' or 'TM').
    - field_string: Field component to extract ('Ex', 'Ey', etc.).
    - pml: Perfectly matched layer settings.
 
    Returns:
    - n_eff: Effective index of the computed mode.
    - field: Field distribution of the computed mode.
    """
 
    try:
        window_width = core_width + 2*clad_width # [nm]
        window_height = core_height + clad_height*2 # [nm]
 
        em.settings(
            wavelength = wavelength,
            x_resolution = dx,
            y_resolution = dy,
            window_width = window_width,
            window_height = window_height,
            num_modes = nr_modes,
            max_effective_index = max_index,
            pml_NSEW_bool = [0,0,0,0],
            boundary_condition = BC,
            background_material = 'SiO2')
 
        em.shape(name='substrate',
                material='SiO2',
                width=window_width,
                height=clad_height)
 
        em.shape(name='ALD_layer',
                material='Al2O3',
                width=window_width,
                height=10)
        em.shape(name='core',
                material='InGaP',
                width=window_width,
                height=core_height,
                etch_depth=core_height-0,
                mask=core_width, sidewall_angle=20)
 
        em.shape(name='ALD_layer3',
                material='Al2O3',
                height=10,
                shape_type='conformal')
 
        em.FDM()
        if plot_bool:
            em.plot()
 
    except Exception as e:
        print('Failed',e)
        traceback.print_exc()
        em.close(save=False)
 
    result = {'n_eff': em.get('effective_index')}
    if return_fields:
        result['field'] = em.get_fields(field_string).field
 
    return result
 
 
def simulate_dual_waveguide(em, 
                            wavelength = 2400, 
                            core_height = 136, 
                            core1_width = 1400, 
                            core2_width = 1400, 
                            gap = 500, 
                            clad_width = 2000, 
                            clad_height = 2000, 
                            dx = 5, 
                            dy = 4, 
                            nr_modes = 2, 
                            max_index = 0, 
                            BC = '00', 
                            field_string = 'Ex', 
                            plot_bool=False, 
                            return_fields=False):
    """
    Simulates an InGaP waveguide using the emodeconnection library.
 
    Parameters:
    - em: An instance of the emodeconnection class.
    - wavelength: Wavelength of light in nm.
    - core_height: Height of the waveguide core in nm.
    - core_widths: Width of the waveguide cores in nm.
    - gap: Gap between the waveguide cores
    - clad_width: Width of the cladding in nm.
    - clad_height: Height of the cladding in nm.
    - dx: Resolution in x-direction in nm.
    - dy: Resolution in y-direction in nm.
    - nr_modes: Number of modes to compute.
    - max_index: Maximum index for mode computation.
    - BC: Boundary condition ('TE' or 'TM').
    - field_string: Field component to extract ('Ex', 'Ey', etc.).
    - pml: Perfectly matched layer settings.
 
    Returns:
    - n_eff: Effective index of the computed mode.
    - field: Field distribution of the computed mode.
    """
 
    try:
        window_width = core1_width + core2_width + gap + 2 * clad_width # [nm]
        window_height = core_height + 2 * clad_height # [nm]
 
        em.settings(
            wavelength = wavelength,
            x_resolution = dx,
            y_resolution = dy,
            window_width = window_width,
            window_height = window_height,
            num_modes = nr_modes,
            max_effective_index = max_index,
            pml_NSEW_bool = [0,0,0,0],
            boundary_condition = BC,
            background_material = 'SiO2')
 
        em.shape(name='substrate',
                material='SiO2',
                width=window_width,
                height=clad_height)
 
        em.shape(name='ALD_layer',
                material='Al2O3',
                width=window_width,
                height=10)
 
        em.shape(name='core',
                material='InGaP',
                width=window_width,
                height=core_height,
                etch_depth=core_height,
                mask=[core1_width, core2_width], mask_offset=[-core1_width/2-gap/2, core2_width/2+gap/2], sidewall_angle=20)
 
        em.shape(name='ALD_layer2',
                material='Al2O3',
                height=10,
                shape_type='conformal')
 
        em.FDM()
        if plot_bool:
            em.plot()
 
    except Exception as e:
        print('Failed',e)
        traceback.print_exc()
        em.close(save=False)
 
    result = {'n_eff': em.get('effective_index')}
    if return_fields:
        result['field'] = em.get_fields(field_string).field
 
    return result

em = emode.EMode() 
result = simulate_waveguide(em, 
                            wavelength=1550,     # nm                            
                            core_height=136,     # nm                           
                            core_width=1200,     # nm                            
                            clad_width=1200,     # nm of SiO2 on each side                            
                            clad_height=1200,    # nm of SiO2 above and below                            
                            dx=10, dy=2,         # nm grid spacing (larger = faster)                            
                            nr_modes=1,                            
                            BC='TE',                            
                            plot_bool=True,      # opens EMode's plot window                            
                            return_fields=False)

print('n_eff:', result['n_eff'])

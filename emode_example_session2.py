import traceback
import emodeconnection as emode
import numpy as np
import scipy.constants as consts
from scipy.integrate import odeint

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
                       field_strings = ['Ex'], 
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
        for field_string in field_strings:
            result[field_string] = em.get_fields(field_string).field

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


####################################################################################################
####################################################################################################
####################################################################################################
#Added in session 2

#Needed for kappa calculation
def calculate_kappa(Ep_X, Ep_Y, Ep_Z, Ei_X, Ei_Y, Ei_Z, Es_X, Es_Y, Es_Z, dx, dy):
    """
    Calculates the coupling coefficient (kappa) based on the provided electric field components.
    NOTE: Needs field components along the crystal axes.

    Parameters:
    - Ex_p, Ey_p, Ez_p: Electric field components of the pump mode.
    - Ex_i, Ey_i, Ez_i: Electric field components of the idler mode.
    - Ex_s, Ey_s, Ez_s: Electric field components of the signal mode.
    - dx: Resolution in x-direction in nm.
    - dy: Resolution in y-direction in nm.

    Returns:
    - kappa: The calculated coupling coefficient.
    """

    term1 = (np.conjugate(Es_X) * np.conjugate(Ei_Y) + np.conjugate(Es_Y) * np.conjugate(Ei_X)) * Ep_Z
    term2 = (np.conjugate(Es_X) * Ep_Y + np.conjugate(Es_Y) * Ep_X) * np.conjugate(Ei_Z)
    term3 = (np.conjugate(Ei_X) * Ep_Y + np.conjugate(Ei_Y) * Ep_X) * np.conjugate(Es_Z)

    # Calculate the numerical integral
    front_factor = 106e-12*consts.c*consts.epsilon_0/2 #106 pm/V is the measured d_14 of InGaP
    integral = np.sum(term1 + term2 + term3) * dx * dy
    kappa = front_factor*integral

    return kappa

def rotate_coordinates(E_x, E_y, E_z):
    E_X = 0
    E_Y = 0
    E_Z = 0
    return E_X, E_Y, E_Z

def core_mask(nx, ny, h, w, clad, x0=0.0):
    """Creates a grid mask for the waveguide core
    1 inside the rectangular core, 0 outside. Shape (ny, nx)."""
    W = 2*clad+w
    H = 2*clad+h
    y0 = clad
    X, Y = np.meshgrid(np.linspace(-W/2, W/2, nx), np.linspace(0, H, ny))
    return ((np.abs(X - x0) <= w/2) & (Y >= y0) & (Y <= y0 + h)).astype(float).T


#Needed for solving the coupled amplitude equations
def db_to_m(a):
    return a * np.log(10) / 10 * 1e2
# https://stackoverflow.com/questions/19910189/scipy-odeint-with-complex-initial-values
def odeintz(func, z0, t, **kwargs):
    """An odeint-like function for complex valued differential equations."""

    # Disallow Jacobian-related arguments.
    _unsupported_odeint_args = ['Dfun', 'col_deriv', 'ml', 'mu']
    bad_args = [arg for arg in kwargs if arg in _unsupported_odeint_args]
    if len(bad_args) > 0:
        raise ValueError("The odeint argument %r is not supported by "
                         "odeintz." % (bad_args[0],))

    # Make sure z0 is a numpy array of type np.complex128.
    z0 = np.array(z0, dtype=np.complex128, ndmin=1)

    def realfunc(x, t, *args):
        z = x.view(np.complex128)
        dzdt = func(z, t, *args)
        # func might return a python list, so convert its return
        # value to an array with type np.complex128, and then return
        # a np.float64 view of that array.
        return np.asarray(dzdt, dtype=np.complex128).view(np.float64)

    result = odeint(realfunc, z0.view(np.float64), t, **kwargs)

    if kwargs.get('full_output', False):
        z = result[0].view(np.complex128)
        infodict = result[1]
        return z, infodict
    else:
        z = result.view(np.complex128)
        return z

def coupled_amplitude_equations(kappa, length, res, a_p, a_i, a_s, Pp_0, Pi_0, Ps_0):
    """
        Solves the coupled amplitude equations for DFG with wavelength 945/1550/2421 over a length L.
    
        Parameters:
        - Kappa: Previously calculated nonlinear overlap.
        - Length: The length (in m) over which to run the simulation.
        - Res: Resolution in the length.
        - a_p, a_i, a_s: Loss values at each wavelength in db/m.
        - Pp_0, Pi_0, Ps_0: Starting power at each wavelength in W.
    
        Returns:
        - P_p, P_i, P_s: Power calculated over the length of the waveguide.
        """
    
    a_p, a_i, a_s = db_to_m(a_p), db_to_m(a_i), db_to_m(a_s)
    z = np.arange(0, length, res)

    'Function that returns Ca computed from an ODE for a k'
    def model(Z, z):
        w_p = 2 * np.pi / (945 * 1e-9)
        w_i = 2 * np.pi / (1550 * 1e-9)
        w_s = 2 * np.pi / (2421 * 1e-9)

        Zb = Z[0]
        Zg = Z[1]
        Zr = Z[2]

        d_alpha = a_p - a_i - a_s
        d_beta = 0
        d_k = d_beta + 1j * d_alpha / 2

        dZpdz = 1j * 2 * w_p * kappa * Zg * Zr * np.exp(-1j * d_k * z)
        dZidz = 1j * 2 * w_i * kappa * Zb * np.conj(Zr) * np.exp(1j * d_k * z) * np.exp(-a_s * z)
        dZsdz = 1j * 2 * w_s * kappa * Zb * np.conj(Zg) * np.exp(1j * d_k * z) * np.exp(-a_i * z)

        return [dZpdz, dZidz, dZsdz]

    Z0 = np.array([np.sqrt(Pp_0) + 0.0j, np.sqrt(Pi_0) + 0.0j, np.sqrt(Ps_0) + 0.0j])
    Z, infodict = odeintz(model, Z0, (z), full_output=True)

    P_p = (Z[:, 0].real ** 2 + Z[:, 0].imag ** 2) * np.exp(-a_p * z)
    P_i = (Z[:, 1].real ** 2 + Z[:, 1].imag ** 2) * np.exp(-a_i * z)
    P_s = (Z[:, 2].real ** 2 + Z[:, 2].imag ** 2) * np.exp(-a_s * z)
    
    return np.array(P_p), np.array(P_i), np.array(P_s)

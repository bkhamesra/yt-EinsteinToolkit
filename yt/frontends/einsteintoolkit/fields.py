"""
EinsteinToolkit-specific fields

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

from yt.fields.magnetic_field import \
    setup_magnetic_field_aliases

density_unit = 'code_mass / code_length**3'

class EinsteinToolkitFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # ADMBase
        ('ADMBASE::alp', ('', ['alpha', 'lapse'], r'$\alpha$')),
        ('ADMBASE::dtalp', ('', ['dtalpha', 'dtlapse', 'time_derivative_lapse'], r'$\partial_t \alpha$')),
        ('ADMBASE::betax', ('', ['shift_x', 'beta_x'], r'$\beta^x$')),
        ('ADMBASE::betay', ('', ['shift_y', 'beta_y'], r'$\beta^y$')),
        ('ADMBASE::betaz', ('', ['shift_z', 'beta_z'], r'$\beta^z$')),
        ('ADMBASE::dtbetax', ('', ['dtbeta_x', 'time_derivative_shift_x'], r'$\partial_t \beta^x$')),
        ('ADMBASE::dtbetay', ('', ['dtbeta_y', 'time_derivative_shift_y'], r'$\partial_t \beta^y$')),
        ('ADMBASE::dtbetaz', ('', ['dtbeta_z', 'time_derivative_shift_z'], r'$\partial_t \beta^z$')),
        ('ADMBASE::gxx', ('', ['metric_xx', 'g_xx'], r'$\gamma_{xx}$')),
        ('ADMBASE::gyy', ('', ['metric_yy', 'g_yy'], r'$\gamma_{yy}$')),
        ('ADMBASE::gzz', ('', ['metric_zz', 'g_zz'], r'$\gamma_{zz}$')),
        ('ADMBASE::gxy', ('', ['metric_xy', 'metric_yx', 'g_xy', 'g_yx'], r'$\gamma_{xy}$')),
        ('ADMBASE::gxz', ('', ['metric_xz', 'metric_xz', 'g_xz', 'g_zx'], r'$\gamma_{xz}$')),
        ('ADMBASE::gyz', ('', ['metric_yz', 'metric_zy', 'g_yz', 'g_zy'], r'$\gamma_{yz}$')),
        ('ADMBASE::kxx', ('', ['extrinsic_curvature_xx', 'k_xx'], r'$K_{xx}$')),
        ('ADMBASE::kyy', ('', ['extrinsic_curvature_yy', 'k_yy'], r'$K_{yy}$')),
        ('ADMBASE::kzz', ('', ['extrinsic_curvature_zz', 'k_zz'], r'$K_{zz}$')),
        ('ADMBASE::kxy', ('', ['extrinsic_curvature_xy', 'k_xy', 'extrinsic_curvature_yx', 'k_yx'], r'$K_{xy}$')),
        ('ADMBASE::kxz', ('', ['extrinsic_curvature_xz', 'k_xz', 'extrinsic_curvature_zx', 'k_zx'], r'$K_{xz}$')),
        ('ADMBASE::kyz', ('', ['extrinsic_curvature_yz', 'k_yz', 'extrinsic_curvature_zy', 'k_zy'], r'$K_{zy}$')),
        # TmunuBase
        ('TMUNUBASE::eTtt', (density_unit, ['T_tt', 'stress_energy_scalar'], r'$T_{tt}$')),
        ('TMUNUBASE::eTtx', (density_unit, ['T_tx', 'T_xt',  'stress_energy_vector_x'], r'$T_{tx}$')),
        ('TMUNUBASE::eTty', (density_unit, ['T_ty', 'T_yt', 'stress_energy_vector_y'], r'$T_{ty}$')),
        ('TMUNUBASE::eTtz', (density_unit, ['T_tz', 'T_zy', 'stress_energy_vector_z'], r'$T_{tz}$')),
        ('TMUNUBASE::eTxx', (density_unit, ['T_xx', 'stress_energy_tensor_xx'], r'$T_{xx}$')),
        ('TMUNUBASE::eTxy', (density_unit, ['T_xy', 'T_yx', 'stress_energy_tensor_xy', 'stress_energy_tensor_yx'], r'$T_{xy}$')),
        ('TMUNUBASE::eTxz', (density_unit, ['T_xz', 'T_zx', 'stress_energy_tensor_xz', 'stress_energy_tensor_zx'], r'$T_{xz}$')),
        ('TMUNUBASE::eTyy', (density_unit, ['T_yy', 'stress_energy_tensor_yy'], r'$T_{yy}$')),
        ('TMUNUBASE::eTyz', (density_unit, ['T_yz', 'T_zy', 'stress_energy_tensor_yz', 'stress_energy_tensor_zy'], r'$T_{yz}$')),
        ('TMUNUBASE::eTzz', (density_unit, ['T_zz', 'stress_energy_tensor_zz'], r'$T_{zz}$')),
        # HydroBase
        ('HYDROBASE::rho', (density_unit, ['density', 'rho'], r'$\rho$')),
        ('HYDROBASE::press', ('code_mass/(code_length*code_time**2)', ['pressure'], r'$P$')),
        ('HYDROBASE::velx', ('c', ['velocity_x'], r'$v^x$')),
        ('HYDROBASE::vely', ('c', ['velocity_y'], r'$v^y$')),
        ('HYDROBASE::velz', ('c', ['velocity_z'], r'$v^z$')),
        ('HYDROBASE::vel[0]', ('c', ['velocity_x'], r'$v^x$')),
        ('HYDROBASE::vel[1]', ('c', ['velocity_y'], r'$v^y$')),
        ('HYDROBASE::vel[2]', ('c', ['velocity_z'], r'$v^z$')),
        ('HYDROBASE::eps', ('c**2', ['eps', 'internal_energy'], r'$\epsilon$')),
        ('HYDROBASE::Bvec[0]', ('', ['magnetic_field_x'], r'$B^x$')),
        ('HYDROBASE::Bvec[1]', ('', ['magnetic_field_y'], r'$B^y$')),
        ('HYDROBASE::Bvec[2]', ('', ['magnetic_field_z'], r'$B^z$')),
        ('HYDROBASE::Avec[0]', ('', ['vector_potential_x', 'magnetic_potential_x'], r'$A^x$')),
        ('HYDROBASE::Avec[1]', ('', ['vector_potential_y', 'magnetic_potential_y'], r'$A^y$')),
        ('HYDROBASE::Avec[2]', ('', ['vector_potential_z', 'magnetic_potential_z'], r'$A^z$')),
        ('HYDROBASE::Aphi', ('', ['scalar_potential', 'electric_potential'], r'$\phi$')),
        # IllinoisGRMHD
        ('ILLINOISGRMHD::rho_star', ('', ['rho_star'], r'$\rho^{*}$')),
        ('ILLINOISGRMHD::tau', ('', ['mhd_tau'], r'$\widetilde{\tau}$')),
        ('ILLINOISGRMHD::mhd_st_x', ('', ['mhd_st_x'], r'$\widetilde{S}_x$')),
        ('ILLINOISGRMHD::mhd_st_y', ('', ['mhd_st_y'], r'$\widetilde{S}_y$')),
        ('ILLINOISGRMHD::mhd_st_z', ('', ['mhd_st_z'], r'$\widetilde{S}_z$')),
        ('ILLINOISGRMHD::rho_b', (density_unit, ['rho_b', 'rho', 'density'], r'$\rho$')),
        ('ILLINOISGRMHD::P', ('code_mass/(code_length*code_time**2)', ['pressure'], r'$P$')),
        ('ILLINOISGRMHD::vx', ('c', ['velocity_x'], r'$v^x$')),
        ('ILLINOISGRMHD::vy', ('c', ['velocity_y'], r'$v^y$')),
        ('ILLINOISGRMHD::vz', ('c', ['velocity_z'], r'$v^z$')),
        ('ILLINOISGRMHD::Bx', ('', ['magnetic_field_x'], r'$B^x$')),
        ('ILLINOISGRMHD::Bxy', ('', ['magnetic_field_y'], r'$B^y$')),
        ('ILLINOISGRMHD::Bxz', ('', ['magnetic_field_z'], r'$B^z$')),
        # WeylScal
        ('WEYLSCAL4::Psi4r', ('', ['psi4r', 'psi4_real'], r'$\Re(\Psi_4)$')),
        ('WEYLSCAL4::Psi4i', ('', ['psi4i', 'psi4_imaginary'], r'$\Im(\Psi_4)$')),
        ('WEYLSCAL4::MagPsi4', ('', ['magpsi4', 'psi4_mag'], r'$\left| \Psi_4 \right|$')),
        # Kranc2BSSNChiMatter
        ('KRANC2BSSNCHIMATTER::chi', ('', ['chi'], r'$\chi$')),
        # CalculateExpansion
        ('CALCULATEEXPANSION::expansion[0]', ('', ['expansion_0'], r'$\Theta_0$')),
        ('CALCULATEEXPANSION::expansion[1]', ('', ['expansion_1'], r'$\Theta_1$')),
        ('CALCULATEEXPANSION::expansion[2]', ('', ['expansion_2'], r'$\Theta_2$')),
        # EMevoPlasmaChi
        ('EMEVOPLASMACHI::Bff1', ('', ['B_x', 'magnetic_field_x'], r'$B^x$')),
        ('EMEVOPLASMACHI::Bff2', ('', ['B_y', 'magnetic_field_y'], r'$B^y$')),
        ('EMEVOPLASMACHI::Bff3', ('', ['B_z', 'magnetic_field_z'], r'$B^z$')),
        ('EMEVOPLASMACHI::Eff1', ('', ['E_x', 'electric_field_x'], r'$E^x$')),
        ('EMEVOPLASMACHI::Eff2', ('', ['E_y', 'electric_field_y'], r'$E^y$')),
        ('EMEVOPLASMACHI::Eff3', ('', ['E_z', 'electric_field_z'], r'$E^z$')),
        ('EMEVOPLASMACHI::charge', ('', ['em_rho', 'charge_density'], r'$\rho$')),
        ('EMEVOPLASMACHI::current1', ('', ['em_j_x', 'current_density_x'], r'$j^x$')),
        ('EMEVOPLASMACHI::current2', ('', ['em_j_y', 'current_density_y'], r'$j^y$')),
        ('EMEVOPLASMACHI::current3', ('', ['em_j_z', 'current_density_z'], r'$j^z$')),
        ('EMEVOPLASMACHI::EdotB', ('', ['electric_dot_magnetic', 'e_dot_b'], r'$E\cdot B$')),
        ('EMEVOPLASMACHI::B2mE2', ('', ['em_b2_minus_e2'], r'$B^2 - E^2$')),
        ('EMEVOPLASMACHI::EMflux', ('', ['em_flux'], r'$\frac{1}{4\pi} \left| \phi_2 \right|^2$')),
    )

    known_particle_fields = ()

    def __init__(self, ds, field_list):
        super(EinsteinToolkitFieldInfo, self).__init__(ds, field_list)

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(self, 'einsteintoolkit', ['EMEVOPLASMACHI::Bff1', 'EMEVOPLASMACHI::Bff2', 'EMEVOPLASMACHI::Bff3'])

    def setup_particle_fields(self, ptype):
        super(EinsteinToolkitFieldInfo, self).setup_particle_fields(ptype)
        # This will get called for every particle type.

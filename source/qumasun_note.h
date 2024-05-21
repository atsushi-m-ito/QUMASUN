#pragma once

namespace QUMASUN {
    static const char* note_def_energy =
        "Note: definitions of energies\n"
        "  Etot = Ekin + E_PP_local_HR + Ehart + Exc + E_PP_nonlocal + Enn(Vlocal*rho_nucl-self_v2)HR + m_Enn_close_correction\n"
        "  E_PP_local_HR = \\int{\\rho V'_local}dV on High-Resolution grid\n"
        "  Eext_nucl_density = \\int{\\rho_nucl V_hart}dV  (non-use)\n"
        "  Eext_nucl_point = \\int{Q_n \\delta(x-R_n) V_hart}dV  (non-use)\n"
        "  E_PP_nonlocal = <psi|V_nonlocal|psi>\n"
        "  Ehart = \\int{\\rho V_hart}dV\n"
        "  Exc = Ex + Ec\n"
        "  Enn(Vlocal*rho_nucl-self_v2)HR = Enn(Vlocal*rho_nucl)HR - Enn(self_v2)HR\n"
        "  Enn(Vlocal*rho_nucl)HR = \\int{\\rho_nucl Vlocal}dV on High-Resolution grid\n"
        "  Enn(self_v2)HR = \\int {\\rho_nucl Vlocal(PP-original)}dV (r < r_cut) on High-Resolution grid\n"
        "    which is self interaction of nuclei estimated as an independent state\n"
        "  Enn(direct) = QQ/r,\n"
        "    where r is distance of two nuclei. (non-use)\n"
        "  m_Enn_close_correction = ZZ/r - \\int{\\rho_nucl Vlocal}dV in ellipsoidal grid,\n"
        "    which is correction of two nuclei when r < core cutoff.\n"
        "  \\rho_nucl is corresponding charge of V_local,\n"
        "    and is converted from radial grid to Cartesian grid.\n"
        "  Vlocal on Cartesian grid is the solution of the Poisson equation from rho_nucl,\n"
        "    and is including the effect of periodic boundary,\n"
        "    while Vlocal(PP-original) is original V_local potential of PP files\n";

}

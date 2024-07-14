title: Appendix
---

[TOC]

## Non-dimensional Navier Stokes equations

To illustrate the roles of various terms in the governing equations, we present here the non-dimensionalized governing compressible Navier Stokes equations with particles. We define the non-dimensional variables as:

\begin{equation}
x_i^* = \frac{x_i}{L_o}, \quad y_i^* = \frac{y_i}{L_o}, \quad t^* = \frac{tu_o}{L_o}, \quad \rho^* = \frac{\rho}{\rho_o}, \quad u_i^* = \frac{u_i}{u_o}
\end{equation}

\begin{equation}
v_i^* = \frac{v_i}{u_o}, \quad p^* = \frac{p}{\rho_o u_o^2}, \quad e^* = \frac{e}{u_o^2}, \quad T^* = \frac{T}{T_o}, \quad T_p^* = \frac{T_p}{T_o}, 
\end{equation}

\begin{equation}
\mu^* = \frac{\mu}{\mu_o}, \quad \kappa^* = \frac{\kappa}{\kappa_o}, \quad g^*_i = \frac{g_i}{g_o}	
\end{equation}

where the subscript \(o\) denotes a reference value. Under these scalings, the Navier Stokes equations become:

\begin{equation}
\label{eq:NS_cont_nd}
\frac{\partial \rho^*}{\partial t^*} + \frac{\partial \rho^* u_j^*}{\partial x_j^*} = 0	
\end{equation}

\begin{equation}
\label{eq:NS_mom_nd}
\begin{split}
\frac{\partial \rho^* u_i^*}{\partial t^*} + \frac{\partial \rho^* u_i^* u_j^*}{\partial x_j^*} = & - \frac{\partial p^*}{\partial x_i^*} + \frac{1}{Re}\frac{\partial \tau_{ij}^*}{\partial x_j^*} + \frac{1}{Fr^2} g_i^* - \\
& \beta\frac{\mu^* \phi_m}{N_p St}  \sum_{n=1}^{N_p} (u_i^*-v_{i,n}^*)\delta(\mathbf{x}^* - \mathbf{y}_n^* )	
\end{split}
\end{equation}

\begin{equation}
\label{eq:NS_en_nd}
\begin{split}
\frac{\partial \rho^* e^*}{\partial t^*} + \frac{\partial u_j^* (\rho^* e^* + p^*)}{\partial x_j^*} = & \frac{1}{Re}\left[\frac{\partial \tau_{ij}^* u_j^*}{\partial x_i^*} + \frac{1}{(\gamma-1)Pr M_o^2}\frac{\partial}{\partial x_j^*}\left( k^* \frac{\partial T^*}{\partial x_j^*}\right)\right]+ \\
& \beta \frac{\phi_m}{3N_p} \frac{Nu}{(\gamma-1)Pr M_o^2St}\sum_{n=1}^{N_p} (T_{p,n}^* - T^*) \delta(\mathbf{x}^* - \mathbf{y}_n^* )
\end{split}
\end{equation}

Where \(Re=\frac{\rho_o L_o u_o}{\mu_o}\) is the Reynolds number, \(Fr=\frac{u_o}{\sqrt{g_o L_o}}\) is the Froude number, \(\phi_m=\frac{m_p N_p}{\rho_o L_o^3}\), \(St=\frac{\tau_p}{\tau_f}\) is the Stokes number which for Stokesian particles is \(St=\frac{\rho_p D_p^2 u_o}{18 L_o \mu_o}\), the Prandtl number (which is assumed constant and equal to 0.72) \(Pr=\frac{\mu_o c_p}{k_o}\), the Nusselt number \(Nu=\frac{h D_p}{k_o}=2\) (for Stokesian particles) and  the Mach number \(M_o=\frac{u_o}{\sqrt{\gamma R T_o}}\).


## Non-dimensional particle equations

The non-dimensional set of equations for the particles reads:

\begin{equation}
\label{eq:part_motion}
\frac{d y_i^*}{dt^*} = v_i^*, \quad \frac{d v_i^*}{dt^*} = \frac{\mu^*}{St}\left(u_i^* - v_i^*\right) + \frac{1}{Fr^2}g_i^*
\end{equation}

\begin{equation}
\label{eq:part_energy}
\frac{d T_p^*}{dt^*} = \frac{\gamma}{3\Gamma St Pr} \left(I_o^* - Nu(T_p^*-T^*)\right)
\end{equation}

where \(I_o^*=\frac{I_o D_p}{4k_oT_o} \) and \(\Gamma=c_{v,p}/c_{v}\) is the ratio of the particle specific heat capacity to the fluid isochoric specific heat capacity.

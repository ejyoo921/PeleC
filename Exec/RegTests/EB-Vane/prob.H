#ifndef PROB_H
#define PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  state(i, j, k, URHO) = prob_parm.rho;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = prob_parm.rho * prob_parm.eint;
  state(i, j, k, UEDEN) = state(i, j, k, UEINT) + 1/state(i, j, k, URHO) *
      (state(i, j, k, UMX)*state(i, j, k, UMX) +  state(i, j, k, UMY)* state(i, j, k, UMY) + state(i, j, k, UMZ)*state(i, j, k, UMZ));
  state(i, j, k, UTEMP) = prob_parm.T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = prob_parm.rho * prob_parm.massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
    // OUTFLOW BCs

    // Initialize all quantities with FOExtrap values
    amrex::Real rho = s_int[URHO];
    amrex::Real u = s_int[UMX] / rho;
    amrex::Real v = s_int[UMY] / rho;
    amrex::Real w = s_int[UMZ] / rho;
    amrex::Real eint = s_int[UEINT] / rho;
    amrex::Real T = s_int[UTEMP];
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = s_int[UFS + n] / rho;
    }

    // Following Blazek p 279, eq. 8.23 - reimplemented from EB-C10
    // NOTE: it is assumed that mass fractions are constant across the
    // boundary Interior state (point d)
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real xd = sgn > 0 ? prob_lo[idir] + 0.5 * dx[idir] : prob_hi[idir] - 0.5 * dx[idir];
    const amrex::Real rho_inter = rho;
    const amrex::Real u_inter = u;
    const amrex::Real v_inter = v;
    const amrex::Real w_inter = w;
    const amrex::Real T_inter = T;
    amrex::Real p_inter = 0.0, cs_inter = 0.0;
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2P(rho_inter, T_inter, massfrac, p_inter);
    eos.RTY2Cs(rho_inter, T_inter, massfrac, cs_inter);

    // Boundary state (point b)
    const amrex::Real xb = sgn > 0 ? prob_lo[idir] : prob_hi[idir];
    const amrex::Real pb = prob_parm.p;
    const amrex::Real rhob =
        s_int[URHO] + (pb - p_inter) / (cs_inter * cs_inter);
    const amrex::Real ub = idir == 0 ? u_inter - sgn * (p_inter - pb) / (rho_inter * cs_inter) : u_inter;
    const amrex::Real vb = idir == 1 ? v_inter - sgn * (p_inter - pb) / (rho_inter * cs_inter) : v_inter;
    const amrex::Real wb = idir == 2 ? w_inter - sgn * (p_inter - pb) / (rho_inter * cs_inter) : w_inter;

    // Ghost state (point a). Linear extrapolation from d and b
    rho = (rhob - rho_inter) / (xb - xd) * (x[idir] - xd) + rho_inter;
    const amrex::Real p = (pb - p_inter) / (xb - xd) * (x[idir] - xd) + p_inter;

    eos.RYP2E(rho, massfrac, p, eint);
    eos.REY2T(rho, eint, massfrac, T);

    u = (ub - u_inter) / (xb - xd) * (x[idir] - xd) + u_inter;
    v = (vb - v_inter) / (xb - xd) * (x[idir] - xd) + v_inter;
    w = (wb - w_inter) / (xb - xd) * (x[idir] - xd) + w_inter;

    // Write to s_ext
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * massfrac[n];
    }
}

void pc_prob_close();

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
    static void problem_modify_ext_sources(
                                           amrex::Real time,
                                           amrex::Real dt,
                                           const amrex::MultiFab& state_old,
                                           const amrex::MultiFab& state_new,
                                           amrex::MultiFab& ext_src,
                                           int ng,
                                           amrex::GeometryData const& geomdata,
                                           ProbParmDevice const& prob_parm)
    {
        /* Notes: ext_src contains sources from velocity forcing coming in
           This function should add to rather than overwrite ext_src.
        */
        auto const& fact =
            dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
        auto const& flags = fact.getMultiEBCellFlagFab();

        auto const& Sns = state_new.const_arrays();
        auto const& Farrs = ext_src.arrays();
        auto const& flagarrs = flags.const_arrays();
        const amrex::IntVect ngs(ng);

        const amrex::Real src_radius2 = prob_parm.source_radius*prob_parm.source_radius;
        const amrex::Real src_strength = prob_parm.source_strength;
        const amrex::Real Tsrc = prob_parm.T;

        amrex::ParallelFor(
                           ext_src, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                               if (!flagarrs[nbx](i, j, k).isCovered()) {
                                   const amrex::RealVect x = pc_cmp_loc({AMREX_D_DECL(i, j, k)}, geomdata);
                                   if (x[0]*x[0] + x[1]*x[1] < src_radius2) {
                                       Farrs[nbx](i, j, k, URHO) += src_strength;
                                       Farrs[nbx](i, j, k, UFS) += src_strength;

                                       auto eos = pele::physics::PhysicsType::eos();
                                       amrex::Real e_int;
                                       amrex::Real massfrac[NUM_SPECIES] = {0.0};
                                       massfrac[0] = 1.0;
                                       eos.TY2E(Tsrc, massfrac, e_int);
                                       Farrs[nbx](i, j, k, UEDEN) += e_int * src_strength;
                                   }
                               }
                           });
        amrex::Gpu::synchronize();
    }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;

#endif

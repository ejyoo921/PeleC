#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
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

// EY: for creating an EB geometry
#include "EB.H"
#include <AMReX_EB2_IF_Rotation.H>
#include "Utilities.H"
#include "Geometry.H"

// EY: 
void parse_params(std::unique_ptr<ProbParmDevice>& prob_parm_device, const amrex::Real* problo, const amrex::Real* probhi);

class EBco2snowInlet
  : public pele::pelec::Geometry::Register<EBco2snowInlet>
{
public:
  static std::string identifier() { return "co2snow-inlet"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];

  const amrex::Real rho_p =
    prob_parm.alpha * std::exp(-(x / prob_parm.sigma) * (x / prob_parm.sigma));
  const amrex::Real rho = prob_parm.rho + rho_p;
  const amrex::Real p = prob_parm.p + rho_p * prob_parm.cs * prob_parm.cs;
  const amrex::Real u_p = prob_parm.cs * rho_p / prob_parm.rho;

  amrex::Real T, eint = 0.0;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.massfrac[n];
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(rho, massfrac, p, eint);
  eos.RYP2T(rho, massfrac, p, T);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u_p;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * u_p * u_p);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * prob_parm.massfrac[n];
  }

  //EY: for EB lathe
  // Get A/As
  // amrex::Real AoverAs = define_nozzle_shape(x, prob_parm.throat_height, prob_parm.nozzle_type) / (0.5*prob_parm.throat_height);
  amrex::Real AoverAs = (x > 0) ?
    x / geomdata.ProbHi(0) * (prob_parm.aoveras_out - 1.0) + 1.0
    : x / geomdata.ProbLo(0) * (prob_parm.aoveras_in - 1.0) + 1.0;
}

// AMREX_GPU_DEVICE
// AMREX_FORCE_INLINE
// void
// bcnormal(
//   const amrex::Real* /*x[AMREX_SPACEDIM]*/,
//   const amrex::Real* /*s_int[NVAR]*/,
//   amrex::Real* /*s_ext[NVAR]*/,
//   const int /*idir*/,
//   const int /*sgn*/,
//   const amrex::Real /*time*/,
//   amrex::GeometryData const& /*geomdata*/,
//   ProbParmDevice const& /*prob_parm*/,
//   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
// {
// }
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
  // idir = 1,2 sgn = pm1 -> all the same as the outlet type 1

  if ((idir == 0) && (sgn == 1)) {
    // Inlet types:
    //    0 - FOExtrap pressure (recommended)
    //    1 - Enforce constant pressure
    //    2 - Characteristic extrapolation of pressure (Blazek Eqn 8.32 - 8.36)

    // Inlet conditions - use interior mass fraction
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = s_int[UFS + n] / s_int[URHO];
    }

    // Apply different inlet types for pressure
    amrex::Real p_in = 0.0;
    if (prob_parm.inlet_type == 0) {
      eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac, p_in);
    } else if (prob_parm.inlet_type == 1) {
      p_in = prob_parm.p_inlet;
    } else {
      amrex::Abort("Invalid inlet boundary type specified");
    }

    // Inlet conditions
    amrex::Real T_in = prob_parm.T_inlet;
    amrex::Real rho_in, eint_in;
    eos.PYT2RE(p_in, massfrac, T_in, rho_in, eint_in);
    amrex::Real Cs_in;
    eos.RTY2Cs(rho_in, T_in, massfrac, Cs_in);
    amrex::Real u_in = prob_parm.M_inlet * Cs_in;
    amrex::Real v_in = 0.0;
    amrex::Real w_in = 0.0;

    // Set boundary cells
    s_ext[URHO] = rho_in;
    s_ext[UMX] = rho_in * u_in;
    s_ext[UMY] = rho_in * v_in;
    s_ext[UMZ] = rho_in * w_in;
    s_ext[UEINT] = rho_in * eint_in;
    s_ext[UEDEN] =
      rho_in * (eint_in + 0.5 * (u_in * u_in + v_in * v_in + w_in * w_in));
    s_ext[UTEMP] = T_in;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho_in * massfrac[n];
    }
  } else if ((idir == 0) && (sgn == -1)) {
    // outlet types:
    //    0 - FOExtrap everything
    //    1 - Characteristic Pressure Outlet (recommended)
    //    2 - Fixed Pressure

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

    // If using a pressure-based outlet, compute that
    if (prob_parm.outlet_type == 0) {
      // Nothing further to be done
    } else if (prob_parm.outlet_type == 1) {
      // Following Blazek p 279, eq. 8.23 - reimplemented from EB-C10
      // NOTE: it is assumed that mass fractions are constant across the
      // boundary Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
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
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = prob_parm.p0;
      const amrex::Real rhob =
        s_int[URHO] + (pb - p_inter) / (cs_inter * cs_inter);
      const amrex::Real ub = u_inter + (p_inter - pb) / (rho_inter * cs_inter);
      const amrex::Real vb = v_inter;
      const amrex::Real wb = w_inter;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_inter) / (xb - xd) * (x[0] - xd) + rho_inter;
      const amrex::Real p = (pb - p_inter) / (xb - xd) * (x[0] - xd) + p_inter;

      eos.RYP2E(rho, massfrac, p, eint);
      eos.REY2T(rho, eint, massfrac, T);

      u = (ub - u_inter) / (xb - xd) * (x[0] - xd) + u_inter;
      v = (vb - v_inter) / (xb - xd) * (x[0] - xd) + v_inter;
      w = (wb - w_inter) / (xb - xd) * (x[0] - xd) + w_inter;
    } else if (prob_parm.outlet_type == 2) {
      amrex::Real p = prob_parm.p0;
      auto eos = pele::physics::PhysicsType::eos();
      eos.PYT2RE(p, massfrac, T, rho, eint);
    } else {
      amrex::Abort("Invalid outlet boundary type specified");
    }

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
}

void pc_prob_close();

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  static amrex::Box the_same_box(const amrex::Box& b) { return b; }

  static void add_problem_derive(
    amrex::DeriveList& derive_lst, const amrex::DescriptorList& desc_lst)
  {
    // Add derives as follows and define the derive function below:
    derive_lst.add(
      "rhoerror", amrex::IndexType::TheCellType(), 1, pc_rhoerror,
      the_same_box);
    derive_lst.addComponent("rhoerror", desc_lst, State_Type, 0, NVAR);
  }

  static void pc_rhoerror(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* /*bcrec*/,
    int /*level*/)
  {
    auto const dat = datfab.array();
    auto arr = derfab.array();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
      geomdata.ProbLoArray();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
      geomdata.CellSizeArray();

    ProbParmDevice const* prob_parm = PeleC::d_prob_parm_device;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const amrex::Real rho_exact =
        prob_parm->rho +
        prob_parm->alpha * std::exp(
                             -((x - prob_parm->cs * time) / prob_parm->sigma) *
                             ((x - prob_parm->cs * time) / prob_parm->sigma));
      arr(i, j, k, 0) = dat(i, j, k, URHO) - rho_exact;
    });
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;

#endif

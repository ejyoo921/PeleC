#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p", PeleC::h_prob_parm_device->p);
    pp.query("T", PeleC::h_prob_parm_device->T);
    pp.query("source_strength", PeleC::h_prob_parm_device->source_strength);
    pp.query("source_radius", PeleC::h_prob_parm_device->source_radius);
  }

  PeleC::h_prob_parm_device->massfrac[0] = 1.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(PeleC::h_prob_parm_device->p,  PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->T, PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->eint);
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}

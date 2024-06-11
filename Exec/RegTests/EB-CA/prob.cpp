#include "prob.H"
// EY: for co2snow inlet
#include <AMReX_EB2_IF_Spline.H>

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
  amrex::ParmParse pp("prob");
  pp.query("rho", PeleC::h_prob_parm_device->rho);
  pp.query("p", PeleC::h_prob_parm_device->p);
  pp.query("alpha", PeleC::h_prob_parm_device->alpha);
  pp.query("sigma", PeleC::h_prob_parm_device->sigma);

  amrex::ParmParse ppeb("eb2");
  ppeb.query("cylinder_radius", PeleC::h_prob_parm_device->radius);

  PeleC::h_prob_parm_device->L = (probhi[0] - problo[0]);

  PeleC::h_prob_parm_device->massfrac[0] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2T(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->massfrac.begin(),
    PeleC::h_prob_parm_device->p, PeleC::h_prob_parm_device->T);
  eos.RPY2Cs(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->p,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->cs);

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rho, p, T, gamma, cs, radius, alpha, sigma"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::h_prob_parm_device->L << "," << PeleC::h_prob_parm_device->rho
    << "," << PeleC::h_prob_parm_device->p << ","
    << PeleC::h_prob_parm_device->T << "," << eos.gamma << ","
    << PeleC::h_prob_parm_device->cs << "," << PeleC::h_prob_parm_device->radius
    << "," << PeleC::h_prob_parm_device->alpha << ","
    << PeleC::h_prob_parm_device->sigma << std::endl;
  ofs.close();
}
}

void
PeleC::problem_post_timestep()
{
  if (verbose <= 0) {
    return;
  }

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real rho_err = 0.0;

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Compute error problem post timestep" << std::endl;
    }

    // Calculate the errors
    for (int lev = 0; lev <= finest_level; lev++) {
      PeleC& pc_lev = getLevel(lev);

      bool local_flag = true;
      rho_err += pc_lev.volWgtSquaredSum("rhoerror", time, local_flag);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    // Get the norm and normalize it
    amrex::Real V = volume.sum(0, false);
    rho_err = std::sqrt(rho_err / V);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME= " << time << " RHO ERROR  = " << rho_err << '\n';

      const int log_index = find_datalog_index("mmslog");
      if (log_index >= 0) {
        std::ostream& data_log2 = parent->DataLog(log_index);

        // Write the quantities at this time
        const int datwidth = 14;
        const int datprecision = 6;
        data_log2 << std::setw(datwidth) << time;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_err;
        data_log2 << std::endl;
      }
    }
  }
}

void
PeleC::problem_post_init()
{
  if (verbose <= 0) {
    return;
  }

  amrex::Real time = state[State_Type].curTime();

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      const int log_index = find_datalog_index("mmslog");
      if (log_index >= 0) {
        std::ostream& data_log2 = parent->DataLog(log_index);
        if (time == 0.0) {
          const int datwidth = 14;
          data_log2 << std::setw(datwidth) << "          time";
          data_log2 << std::setw(datwidth) << "       rho_err";
          data_log2 << std::endl;
        }
      }
    }
  }
}

void
PeleC::problem_post_restart()
{
}

// EY: Parm_parse
void
parse_params(ProbParmDevice* prob_parm_device, const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  amrex::ParmParse ppe("eb2");
  std::string geom_type;
  ppe.get("geom_type", geom_type);

  if (geom_type == "co2snow-inlet") {
    amrex::ParmParse ppn("nozzle");
    amrex::Vector<amrex::Real> yvals;
    amrex::Vector<amrex::Real> xvals;
    ppn.getarr("xvals",xvals);
    ppn.getarr("yvals",yvals);
    amrex::Real offset;
    amrex::Real factor;
    ppn.get("offset",offset);
    ppn.get("factor",factor);
    const int nvals = xvals.size();
    for (int i =0; i < nvals; ++i) {
      yvals[i] = offset + factor * yvals[i];
    }
    int i;
    for (i =0; i < nvals; ++i) {
      if (problo[0] < xvals[i]) {
        break;
      }
    }
    amrex::Real r_in = yvals[i];
    for (i = i; i < nvals-1; ++i) {
      if (yvals[i+1] > yvals[i]) {
        break;
      }
    }
    amrex::Real r_t = yvals[i];
    for (i = i; i < nvals; ++i) {
      if (probhi[0] < xvals[i]) {
        break;
      }
    }
    amrex::Real r_out = yvals[i-1];

    bool do_lathe;
    ppn.get("do_lathe", do_lathe);
    prob_parm_device->aoveras_in = (r_in/r_t);
    prob_parm_device->aoveras_out = (r_out/r_t);
    if (do_lathe) {
      prob_parm_device->aoveras_in *= prob_parm_device->aoveras_in;
      prob_parm_device->aoveras_out *= prob_parm_device->aoveras_out;
    }
    prob_parm_device->geom_type = 2;

  }
  else {
    amrex::Abort("Invalid geom_type specified");
  }
}

void
EBco2snowInlet::build( // We don't do spline here; only lathe
  const amrex::Geometry& geom, const int max_coarsening_level)
{

  amrex::EB2::SplineIF nozzle_upper_surface;
  std::vector<amrex::RealVect> lnpts;

  amrex::ParmParse pp("nozzle");
  amrex::Real offset = 1.0;
  amrex::Real factor = 1.0;
  bool do_spline = false;
  bool do_lathe = true;
  pp.query("offset", offset);
  pp.query("factor", factor);
  pp.query("do_spline", do_spline);
  pp.query("do_lathe", do_lathe);

  amrex::Vector<amrex::Real> xvals;
  amrex::Vector<amrex::Real> yvals;
  pp.getarr("xvals",xvals);
  pp.getarr("yvals",yvals);
  int nvals = xvals.size();
  AMREX_ALWAYS_ASSERT(yvals.size() == nvals);
  for (int i = 0; i < nvals; ++i) {
    yvals[i] = offset + factor * yvals[i];
  }

  for (int ii =0;ii<nvals; ++ii) {
    lnpts.push_back(amrex::RealVect(AMREX_D_DECL(xvals[ii], yvals[ii], 0.0)));
  }
  nozzle_upper_surface.addLineElement(lnpts);

  // lathe only supported in 3D
  AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM==3);
  // revolving around x-axis
  auto nozzle_rotate = amrex::EB2::rotate(nozzle_upper_surface, -0.5*constants::PI(), 2);
  auto nozzle_lathe = amrex::EB2::lathe(nozzle_rotate);
  auto nozzle = amrex::EB2::rotate(nozzle_lathe, -0.5*constants::PI(), 1);

  auto gshop = amrex::EB2::makeShop(nozzle);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);

}
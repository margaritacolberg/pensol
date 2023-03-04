// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// main.cc simulates the dynamics of protein folding using a coarse-grained
// model; each amino acid in the protein is represented by a bead, and the
// beads are connected by local and nonlocal bonds; the dynamics is
// event-driven, and a bath of solvent particles facilitates bond forming and
// breaking events

#include "hardspheres.h"
#include "json.hpp"
#include "snapshot.h"
#include "solvent.h"
#include "writer.h"
#include <H5Cpp.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace po = boost::program_options;

void initialize_pos(System &sys, Random &mt, const Param &p, const Box &box,
                    UpdateConfig &update_config,
                    std::optional<std::string> input_name,
                    std::optional<std::string> snapshot_name) {
  if (snapshot_name && boost::filesystem::exists(*snapshot_name)) {
    // overwrite existing entries in pos vectors with read-in
    // values from hdf5 file
    read_snapshot(*snapshot_name, sys.pos, mt, update_config);
  } else if (input_name) {
    read_input(*input_name, sys.pos);
    init_update_config(sys.pos, update_config, box, p.rc2, p.transient_bonds);
  } else {
    init_pos(sys.pos, box, mt, p);
  }
}

void initialize_system(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, Cells &cells_bead,
                       EventQueue &event_queue) {
  if (!check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                        p.nnear_max2)) {
    throw std::runtime_error("local beads overlap");
  }

  if (!check_nonlocal_dist(sys.pos, box, p.rc2, p.rh2, p.permanent_bonds)) {
    throw std::runtime_error("nonlocal beads overlap");
  }

  if (!(cells_bead.ncell >= 4)) {
    throw std::invalid_argument("bead ncell must be at least 4");
  }

  if (!(cells_bead.lcell >= p.rc)) {
    throw std::invalid_argument("bead lcell must be at least rc");
  }

  init_vel(sys.vel, mt, p.temp, p.m);

  // split the box into smaller cells, and store the
  // beads in each cell
  init_cells(sys.pos, box, cells_bead);

  // fill priority queue with cell crossings of all particles
  init_cell_events(sys.pos, sys.vel, p.nbeads, box, sys.counter, event_queue,
                   sys.times, cells_bead);

  // fill priority queue with nearest bond events between all particle
  // pairs
  init_nearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                           event_queue, sys.times, p.near_min2, p.near_max2);

  // fill priority queue with next-nearest bond events between all
  // particle pairs
  init_nnearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                            event_queue, sys.times, p.nnear_min2, p.nnear_max2);

  // fill priority queue with collisions of all particle pairs (executed
  // once, with an initial prediction that fills the entire priority
  // queue)
  add_events_for_all_beads(sys.pos, sys.vel, p.nbeads, p.rh2, p.rc2, box,
                           sys.counter, event_queue, sys.times, cells_bead,
                           p.transient_bonds, p.permanent_bonds, update_config,
                           p.max_nbonds);
}

void run_step(System &sys, const Param &p, const Box &box,
              UpdateConfig &update_config,
              UpdateConfigWriter &update_config_writer, uint64_t countevents,
              uint64_t validevents, double wall_time, Cells &cells_bead,
              EventQueue &event_queue, double step_time) {
  // while events are occurring in step_time,
  while (!event_queue.empty()) {
    // access the minimum time to the next collision, and the
    // indices and collision counters of the beads
    // associated with this collision
    const Event event = event_queue.top();

    if (std::visit(
            [=](auto &&ev) {
              // check for monotonically increasing event times
              assert(ev.t >= wall_time);
              return ev.t > step_time;
            },
            event))
      break;

    countevents++;
    event_queue.pop();

    // process collision or cell crossing event
    if (std::visit(
            [&](auto &&ev) {
              wall_time = ev.t;
              LOG_DEBUG("wall time " << wall_time);
              return process_event(ev, sys, p, box, event_queue, cells_bead,
                                   update_config, update_config_writer);
            },
            event))
      validevents++;
  }

  // update positions at the moment the last collision in p.del_t
  // occurred
  for (unsigned int i = 0; i < p.nbeads; i++) {
    update_pos(sys.pos[i], sys.vel[i], sys.times[i], step_time);
    assert(check_overlap(i, sys.pos, sys.vel, sys.times, p.sigma_bb2, box));
  }

  // check bond distances of nearest and next-nearest beads
  assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2));

  // check bond distances of nonlocal beads
  assert(check_nonlocal_dist(sys.pos, box, p.rc2, p.rh2, p.permanent_bonds));
}

// equilibrate initial configuration, before carrying out bond forming or
// breaking events
void run_trajectory_eq(System &sys, const Param &p, const Box &box,
                       UpdateConfig &update_config,
                       UpdateConfigWriter &update_config_writer,
                       Cells &cells_bead, EventQueue &event_queue,
                       uint64_t countevents, uint64_t validevents, Random &mt,
                       double wall_time) {
  // assume that the entire time during which the beads are undergoing events
  // can be divided into intervals called p.del_t; the total number of such
  // intervals is p.nsteps (thus, the variable called step marks the
  // intervals until the p.nsteps-th interval is reached); each iteration of
  // the BIIIIIIIG loop will run until update_time and then dump the output
  // to the hdf5 file
  for (unsigned int step = 0; step < p.nsteps; step++) {
    std::cout << "step = " << step << std::endl;

    // the current time interval the events are occurring in
    double step_time = step * p.del_t;

    run_step(sys, p, box, update_config, update_config_writer, countevents,
             validevents, wall_time, cells_bead, event_queue, step_time);

    unsigned int freq_coll = p.del_t_coll / p.del_t;
    if (step % freq_coll == 0) {
      collide_sol(sys, mt, p);

      // update counters of the colliding particles
      for (unsigned int i = 0; i < p.nbeads; i++) {
        sys.counter[i]++;
      }

      init_cell_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                       event_queue, sys.times, cells_bead);
      init_nearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                               event_queue, sys.times, p.near_min2,
                               p.near_max2);
      init_nnearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                                event_queue, sys.times, p.nnear_min2,
                                p.nnear_max2);
      add_events_for_all_beads(sys.pos, sys.vel, p.nbeads, p.rh2, p.rc2, box,
                               sys.counter, event_queue, sys.times, cells_bead,
                               p.transient_bonds, p.permanent_bonds,
                               update_config, p.max_nbonds);

      assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2,
                              p.nnear_min2, p.nnear_max2));
      assert(
          check_nonlocal_dist(sys.pos, box, p.rc2, p.rh2, p.permanent_bonds));
    }

    // update time
    wall_time = step_time;
  }
}

void run_trajectory(System &sys, const Param &p, const Box &box,
                    UpdateConfig &update_config,
                    UpdateConfigWriter &update_config_writer, Cells &cells_bead,
                    EventQueue &event_queue, PosWriter &pos_writer,
                    VelWriter &vel_writer, ConfigWriter &config_writer,
                    std::set<Config> &store_config, ConfigInt &store_config_int,
                    std::vector<uint64_t> &config_sum, uint64_t countevents,
                    uint64_t validevents, Random &mt, double wall_time,
                    unsigned int iter, const H5::DataSpace &mem_space,
                    H5::DataSpace &file_space, H5::DataSet &dataset_pos) {
  for (unsigned int step = 0; step < p.nsteps; step++) {
    std::cout << "step = " << step << std::endl;

    double step_time = step * p.del_t;

    run_step(sys, p, box, update_config, update_config_writer, countevents,
             validevents, wall_time, cells_bead, event_queue, step_time);

    config_sum[step] += update_config.config;

    // output positions
    const hsize_t count[3] = {1, p.nbeads, 3};
    const hsize_t start[3] = {step + (p.nsteps * iter), 0, 0};
    file_space.selectHyperslab(H5S_SELECT_SET, count, start);
    dataset_pos.write(&sys.pos[0], H5::PredType::NATIVE_DOUBLE, mem_space,
                      file_space);

    if (step % p.write_step == 0) {
      // store the integer of the configuration and the time
      // of the event
      store_config_int.emplace_back(update_config.config);
      update_config_writer.config_int.emplace_back(update_config.config);
      update_config_writer.config_time.emplace_back(step_time);
      vel_writer.append(sys.vel);

      // if a transient bond forms, check if configuration has
      // been previously visited by comparing configuration to
      // set of saved configurations; if not visited, save
      // configuration to the set and write positions of beads
      // to file
      if ((update_config.config == 1 &&
           store_config.insert(update_config.config).second) ||
          p.transient_bonds.get_nbonds() == 0) {
        pos_writer.append(sys.pos);
        config_writer.append(update_config.config);
      }
    }

    update_config_writer.append();
    update_config_writer.clear();

    unsigned int freq_coll = p.del_t_coll / p.del_t;
    if (step % freq_coll == 0) {
      collide_sol(sys, mt, p);

      for (unsigned int i = 0; i < p.nbeads; i++) {
        sys.counter[i]++;
      }

      init_cell_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                       event_queue, sys.times, cells_bead);
      init_nearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                               event_queue, sys.times, p.near_min2,
                               p.near_max2);
      init_nnearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                                event_queue, sys.times, p.nnear_min2,
                                p.nnear_max2);
      add_events_for_all_beads(sys.pos, sys.vel, p.nbeads, p.rh2, p.rc2, box,
                               sys.counter, event_queue, sys.times, cells_bead,
                               p.transient_bonds, p.permanent_bonds,
                               update_config, p.max_nbonds);

      assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2,
                              p.nnear_min2, p.nnear_max2));
      assert(
          check_nonlocal_dist(sys.pos, box, p.rc2, p.rh2, p.permanent_bonds));
    }

    wall_time = step_time;
  }
}

int main(int argc, char *argv[]) {
  // restart program from point of interruption
  // see Prus, 2002, Boost.Program_options doc
  //
  // declare the supported options
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")(
      "json-file", po::value<std::string>()->required(), "json")(
      "output-file", po::value<std::string>()->required(),
      "hdf5 output")("input-file", po::value<std::string>(), "hdf5 input")(
      "snapshot-file", po::value<std::string>(), "hdf5 snapshot");

  po::positional_options_description pod;
  pod.add("json-file", 1);
  pod.add("output-file", 1);

  po::variables_map vm;
  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pod).run(),
        vm);

    if (vm.count("help")) {
      std::cout << "Usage: hardspheres [options]... json-file output-file"
                << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 0;
    }

    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "hardspheres: " << e.what() << std::endl;
    std::cerr << "try hardspheres --help" << std::endl;
    return 1;
  }

  std::optional<std::string> input_name;
  if (vm.count("input-file")) {
    input_name = vm["input-file"].as<std::string>();
  }

  std::optional<std::string> snapshot_name;
  if (vm.count("snapshot-file")) {
    snapshot_name = vm["snapshot-file"].as<std::string>();
  }

  const std::string json_name = vm["json-file"].as<std::string>();
  const std::string output_name = vm["output-file"].as<std::string>();

  std::cout << "git commit " << VERSION << std::endl;

  std::ifstream input(json_name);
  nlohmann::json json;
  input >> json;

  const Param p = json;
  Param p_eq = p;
  const Box box{p.length};
  ConfigInt store_config_int;

  std::vector<uint64_t> config_sum(p.nsteps, 0);
  std::vector<double> config_prob(p.nsteps);

  System sys(p.nbeads);

  std::seed_seq seq(p.seeds.begin(), p.seeds.end());
  Random mt(seq);

  // open output file
  H5::H5File file(output_name + ".tmp", H5F_ACC_TRUNC);

  UpdateConfigWriter update_config_writer(file);

  H5::Group group{file.createGroup("unique_config")};
  PosWriter pos_writer{group, "pos", p.nbeads};
  VelWriter vel_writer{group, "vel", p.nbeads};
  ConfigWriter config_writer{group, "config"};
  ConfigProbWriter config_prob_writer{group, "config_prob", p.nsteps};

  std::set<Config> store_config;

  H5::Attribute attr_sigma_bb =
      file.createAttribute("sigma_bb", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
  attr_sigma_bb.write(H5::PredType::NATIVE_DOUBLE, &p.sigma_bb);
  H5::Attribute attr_length =
      file.createAttribute("length", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
  attr_length.write(H5::PredType::NATIVE_DOUBLE, &box.l);
  H5::Attribute attr_nsteps =
      file.createAttribute("nsteps", H5::PredType::NATIVE_UINT, H5S_SCALAR);
  attr_nsteps.write(H5::PredType::NATIVE_UINT, &p.nsteps);

  p.nonlocal_bonds.write_hdf5(file, "nonlocal_bonds");
  p.transient_bonds.write_hdf5(file, "transient_bonds");
  p.permanent_bonds.write_hdf5(file, "permanent_bonds");

  const hsize_t mem_dims[2] = {p.nbeads, 3};
  const hsize_t file_dims[3] = {p.nsteps * p.total_iter, p.nbeads, 3};
  H5::DataSpace mem_space(2, mem_dims);
  H5::DataSpace file_space(3, file_dims);
  H5::DataSet dataset_pos =
      file.createDataSet("pos", H5::PredType::NATIVE_FLOAT, file_space);

  uint64_t countevents = 0;
  uint64_t validevents = 0;

  // the BIIIIIIIG loop
  for (unsigned int iter = 0; iter < p.total_iter; iter++) {
    // reset bead clocks, counters, and wall time
    for (unsigned int i = 0; i < p.nbeads; i++) {
      sys.times[i] = 0.0;
      sys.counter[i] = 0.0;
    }

    double wall_time = 0.0;

    Cells cells_bead{p.ncell, p.length / p.ncell};

    UpdateConfig update_config_eq;
    UpdateConfig update_config;

    EventQueue event_queue_eq;
    EventQueue event_queue;

    // reset configuration count
    store_config_int.clear();

    p_eq.transient_bonds = p.transient_bonds_eq;
    p_eq.nsteps = p.nsteps_eq;

    initialize_pos(sys, mt, p_eq, box, update_config_eq, input_name,
                   snapshot_name);
    initialize_system(sys, mt, p_eq, box, update_config_eq, cells_bead,
                      event_queue_eq);

    run_trajectory_eq(sys, p_eq, box, update_config_eq, update_config_writer,
                      cells_bead, event_queue_eq, countevents, validevents, mt,
                      wall_time);

    for (unsigned int i = 0; i < p.nbeads; i++) {
      sys.times[i] = 0.0;
      sys.counter[i] = 0.0;
    }

    wall_time = 0.0;

    init_update_config(sys.pos, update_config, box, p.rc2, p.transient_bonds);
    initialize_system(sys, mt, p, box, update_config, cells_bead, event_queue);

    run_trajectory(sys, p, box, update_config, update_config_writer, cells_bead,
                   event_queue, pos_writer, vel_writer, config_writer,
                   store_config, store_config_int, config_sum, countevents,
                   validevents, mt, wall_time, iter, mem_space, file_space,
                   dataset_pos);

    if (snapshot_name) {
      // write new snapshot to temporary file
      const std::string snapshot_tmp = *snapshot_name + ".tmp";
      write_snapshot(snapshot_tmp, sys.pos, mt, update_config);
      // atomically overwrite old snapshot with new snapshot
      boost::filesystem::rename(snapshot_tmp, *snapshot_name);
    }
  }

  for (unsigned int i = 0; i < p.nsteps; i++) {
    config_prob[i] = double(config_sum[i]) / double(p.total_iter);
  }
  config_prob_writer.append(config_prob);

  H5::Attribute attr_countevents = file.createAttribute(
      "countevents", H5::PredType::NATIVE_UINT64, H5S_SCALAR);
  attr_countevents.write(H5::PredType::NATIVE_UINT64, &countevents);
  H5::Attribute attr_validevents = file.createAttribute(
      "validevents", H5::PredType::NATIVE_UINT64, H5S_SCALAR);
  attr_validevents.write(H5::PredType::NATIVE_UINT64, &validevents);

  file.close();

  std::filesystem::rename(output_name + ".tmp", output_name);

  return 0;
}

void from_json(const nlohmann::json &json, Param &p) {
  p.m = json["m"];
  p.m_sol = json["m_sol"];
  p.sigma_bb = json["sigma_bb"];
  p.sigma_bb2 = p.sigma_bb * p.sigma_bb;
  p.sigma_sb = json["sigma_sb"];
  p.sigma_sb2 = p.sigma_sb * p.sigma_sb;
  p.near_min = json["near_min"];
  p.near_max = json["near_max"];
  p.near_min2 = p.near_min * p.near_min;
  p.near_max2 = p.near_max * p.near_max;
  p.nnear_min = json["nnear_min"];
  p.nnear_max = json["nnear_max"];
  p.nnear_min2 = p.nnear_min * p.nnear_min;
  p.nnear_max2 = p.nnear_max * p.nnear_max;
  p.rh = json["rh"];
  p.rc = json["rc"];
  p.rh2 = p.rh * p.rh;
  p.rc2 = p.rc * p.rc;
  p.eps = json["eps"];
  p.nonlocal_bonds = json["nonlocal_bonds"];
  p.transient_bonds = json["transient_bonds"];
  p.transient_bonds_eq = json["transient_bonds_eq"];
  p.permanent_bonds = json["permanent_bonds"];
  p.tries = json["tries"];
  p.nbeads = json["nbeads"];
  p.length = json["length"];
  p.ncell = json["ncell"];
  p.nsteps = json["nsteps"];
  p.nsteps_eq = json["nsteps_eq"];
  p.total_iter = json["total_iter"];
  p.del_t = json["del_t"];
  p.del_t_coll = json["del_t_coll"];
  p.write_step = json["write_step"];
  p.seeds = json["seeds"].get<std::vector<unsigned int>>();
  p.k_B = json["k_B"];
  p.temp = json["temp"];
  p.nsol = json["nsol"];
  p.max_nbonds = json["max_nbonds"];
  p.angle = json["angle"];
}

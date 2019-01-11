#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <zconf.h>
#include <memory>
#include <map>
#include <cassert>
#include <algorithm>

#include "../include/rules.hpp"
#include "../include/zoltan_fn.hpp"
#include "../include/unloading_model.hpp"
#include "../include/tca_io.hpp"
#include "../include/gif.h"
//#include "../include/geometric_load_balancer.hpp"
#include "zupply.hpp"

#define SPATIAL_DISCRETISATION_X = 7.5 //the average length a conventional vehicle occupies in a closely packed jam (and as such, its width is neglected),
#define TEMPORAL_DISCRETISATION = 1.0 //typical driver’s reaction time

#define START_TIMING(v) \
auto v = MPI_Wtime()

#define RESTART_TIMING(v) \
v = MPI_Wtime() - v

#define STOP_TIMING(v) \
v = MPI_Wtime() - v

#define PAR_START_TIMING(v, comm) \
MPI_Barrier(comm); \
auto v = MPI_Wtime()

#define PAR_RESTART_TIMING(v, comm) \
std::cout <<"Intermediate result " << v << std::endl;\
MPI_Barrier(comm); \
v = MPI_Wtime() - v

#define PAR_STOP_TIMING(v, comm) \
MPI_Barrier(comm); \
v = MPI_Wtime() - v

using namespace std;
using namespace tca;

///TODO: Colormap as function of waiting time
///Parse data afterward

std::tuple<std::vector<Vehicle>, vector<vector<Vehicle *> > > update(const int msx, const int msy,
                                                                     const vector<vector<CA_Cell> > &ca_matrix,
                                                                     const std::vector<Vehicle> &vehicles,
                                                                     const vector<vector<Vehicle *> > &vehicles_map) {
    std::vector<Vehicle> vehicles_new;
    vector<vector<Vehicle *> > vehicles_map_new(msy, std::vector<Vehicle *>(msx, nullptr));
    for (const Vehicle &v : vehicles) {
        deprecated::apply_rule184(msx, msy, ca_matrix, v, vehicles_map, vehicles_new, vehicles_map_new);
    }

    for(auto& v : vehicles_new) {
        int x,y; std::tie(x,y) = v.position;
        std::cout << v << std::endl;
        vehicles_map_new[y][x] = &v;
    }

    return std::make_tuple(vehicles_new, vehicles_map_new);
}

unordered_map<long long, Vehicle> update(const int msx, const int msy,
                                         const unordered_map<long long,  CA_Cell> &ca_matrix,
                                         const unordered_map<long long,  Vehicle> &vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new;
    for (const auto &v : vehicles_map) {
        sequential::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, vehicles_map_new);
    }

    return vehicles_map_new;
}

unordered_map<long long, Vehicle> parallel_update_new_model(const int msx, const int msy,
                                                  const unordered_map<long long,  CA_Cell> &ca_matrix,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map_top,
                                                  const unordered_map<long long,  Vehicle> &remote_vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new, remote_vehicles_map_merged;
    remote_vehicles_map_merged.insert(remote_vehicles_map.begin(), remote_vehicles_map.end());
    remote_vehicles_map_merged.insert(vehicles_map_top.begin(),    vehicles_map_top.end());
    for (const auto &v : vehicles_map) {
        parallel::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, remote_vehicles_map_merged, &vehicles_map_new);
    }

    return vehicles_map_new;
}

unordered_map<long long, Vehicle> parallel_update(const int msx, const int msy,
                                                  const unordered_map<long long,  CA_Cell> &ca_matrix,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map,
                                                  const unordered_map<long long,  Vehicle> &remote_vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new;
    for (const auto &v : vehicles_map) {
        parallel::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, remote_vehicles_map, &vehicles_map_new);
    }

    return vehicles_map_new;
}

void generate_vehicles(const int step, const int my_rank,
                       const int msx, const int msy, Zoltan_Struct* zz,
                       const unordered_map<long long,  CA_Cell> &ca_matrix,
                             unordered_map<long long,  Vehicle> *_vehicles_map) {
    auto& vehicles_map = *_vehicles_map;
    std::vector<std::pair<long long, CA_Cell>> cells;
    CA_Cell c;
    for(long long i = 0; i < msx; i+=1) {
        c = ca_matrix.at(i + 0*msx);
        if(c.source) cells.emplace_back(i + 0*msx, c);
        c = ca_matrix.at(i + (msy-1)*msx);
        if(c.source) cells.emplace_back(i + (msy-1)*msx, c);

    }
    for(long long i = 0; i < msy; i+=1) {
        c = ca_matrix.at(i * msx);
        if (c.source) cells.emplace_back(i * msx, c);
        c = ca_matrix.at((msx - 1) + i * msx);
        if (c.source) cells.emplace_back((msx - 1) + i * msx, c);
    }
    //std::copy_if(ca_matrix.cbegin(), ca_matrix.cend(), std::back_inserter(cells), [](auto cell){return cell.second.source;});
    //add new vehicles
    std::array<double, 2> pos = {0, 0};
    int PE;
    for (auto &cell : cells) {
        std::tie(pos[0], pos[1]) = cell_to_position(msx, msy, cell.first);
        Zoltan_LB_Point_Assign(zz, &pos.front(), &PE);
        if(my_rank == PE && !exists(vehicles_map, cell.first) && cell.second.has_to_generate(step)) {
            vehicles_map[cell.first] = Vehicle(1, vehicles_map.size() + 1, (int) pos[0], (int) pos[1], 1);
        }
    }
}

int main(int argc, char **argv) {
    int wsize, rank;
    srand(1);
    MPI_Init(&argc, &argv);
    MPI_Comm bottom = MPI_COMM_WORLD;
    MPI_Comm_size(bottom, &wsize);
    MPI_Comm_rank(bottom, &rank);
    int SIZE_X, SIZE_Y;
    std::ofstream out;

    unordered_map<long long, CA_Cell> ca_matrix;
    const std::string prefix_fname(argv[1]);
    const int MAX_STEP = std::atoi(argv[2]);
    std::vector<int> road_y_pos;


    SIZE_X = std::atoi(argv[3]);
    SIZE_Y = std::atoi(argv[4]);
    std::tie(ca_matrix, std::ignore, road_y_pos) = generate_random_manhattan(SIZE_X, SIZE_Y);

    //    std::tie(SIZE_X, SIZE_Y) = read_roadfile(argv[3], &ca_matrix);

    create_random_left_sources(3, SIZE_X, SIZE_Y, road_y_pos, &ca_matrix);

    if(!rank) std::cout << "End of map generation\nStarting computations..." << std::endl;

    //

    //vector<vector<CA_Cell> > ca_matrix;

    //vector<Vehicle> vehicles;

    unordered_map<long long, Vehicle> vehicle_matrix;
    //vector<vector<Vehicle *> > vehicle_matrix(SIZE_Y, vector<Vehicle *>(SIZE_X, nullptr));

    auto datatype = Vehicle::register_datatype();
    if(!rank)
        randomize_cars_position(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
    auto vehicles = to_vec(vehicle_matrix);

    if(!rank) std::cout << vehicles.size() << std::endl;
    //std::for_each(vehicles.begin(), vehicles.end(), [](auto v){ std::cout << v << std::endl; });

    //fprint(out, SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
    if(!rank){
        auto img = zzframe( SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
        img.save("out0.jpg");
    }
    auto zz = zoltan_create_wrapper(ENABLE_AUTOMATIC_MIGRATION, MPI_COMM_WORLD);

    //auto remote_data = zoltan_exchange_data(vehicles, zz, datatype.elements_datatype, bottom, recv, sent);
    //std::for_each(remote_data.begin(), remote_data.end(), [](auto v){std::cout << v << std::endl;});

    zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);

    double slope = rank == 0 ? 0.3 : 0.0;
    std::vector<int> incr_cpu;

    auto err = load_balancing::esoteric::get_communicator(slope, rank, bottom, &incr_cpu);

    //if(err) std::cout << err << " no esoteric call" << std::endl;

    std::vector<Vehicle> top_vehicles;
    auto zztop = load_balancing::esoteric::divide_data_into_top_bottom(&vehicles , &top_vehicles, incr_cpu, datatype.elements_datatype, bottom);

    int step = 0;

    while (step < MAX_STEP) {


        MPI_Barrier(bottom);
        PAR_START_TIMING(step_time, bottom);
        int recv, sent;
        /*************************************Start parallel exchange********************************************/
        PAR_START_TIMING(comm_time, bottom);
        auto remote_data = load_balancing::esoteric::exchange(zz, zztop, &vehicles, &top_vehicles, &recv, &sent, incr_cpu, datatype.elements_datatype, bottom, 1.0);
        //auto remote_data =  zoltan_exchange_data(zz, &vehicles, &recv, &sent, datatype.elements_datatype, bottom,  1.2);
        PAR_STOP_TIMING(comm_time, bottom);
        // Stop parallel exchange

        /*************************************Start parallel computation*****************************************/
        PAR_START_TIMING(computation_time, bottom);
        auto vehicle_matrix_remote = to_map(SIZE_X, SIZE_Y, remote_data);
        auto vehicle_matrix_bottom = to_map(SIZE_X, SIZE_Y, vehicles);
        auto vehicle_matrix_top    = to_map(SIZE_X, SIZE_Y, top_vehicles);
        auto updated_vehicle_matrix1 = parallel_update(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_bottom, vehicle_matrix_remote);
        auto updated_vehicle_matrix2 = parallel_update_new_model(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_top, vehicle_matrix_bottom, vehicle_matrix_remote);
        generate_vehicles(step, rank, SIZE_X, SIZE_Y, zz, ca_matrix, &updated_vehicle_matrix1);
        vehicles     = to_vec(updated_vehicle_matrix1);
        top_vehicles = to_vec(updated_vehicle_matrix2);
        PAR_STOP_TIMING(computation_time, bottom);

        // Stop parallel computation

        /********************************Start load balancing and migration**************************************/
        PAR_START_TIMING(migrate_time, bottom);
        load_balancing::esoteric::migrate(zz, zztop, &vehicles, &top_vehicles, incr_cpu, datatype.elements_datatype, bottom );
        //zoltan_migrate_particles(zz, &vehicles, datatype.elements_datatype, bottom);
        //zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);
        PAR_STOP_TIMING(migrate_time, bottom);



        // Stop load balancing and migration
        MPI_Barrier(bottom);
        PAR_STOP_TIMING(step_time, bottom);

        /****************************************Start printing**************************************************/
        if(!rank)
            std::cout << "Time for step " << step
                      << "; [TOT " << step_time
                      << ", CPT " << computation_time
                      << ", COM EXCHANGE " << comm_time
                      << ", COM * MIGRATE " << migrate_time
                      << "] => " << (100*computation_time/step_time)<<"% CPT "
                      << (100*comm_time/step_time)<<"% COM" << std::endl;

        std::vector<Vehicle> all_vehicles;

        if(wsize > 1) {
            gather_elements_on(vehicles, 0, &all_vehicles, datatype.elements_datatype, bottom);
            gather_elements_on(top_vehicles, 0, &all_vehicles, datatype.elements_datatype, bottom);
        }
        else all_vehicles = vehicles;

        if(!rank) {
            auto vehicle_matrix_print = to_map(SIZE_X, SIZE_Y, all_vehicles);
            auto img = zzframe( SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_print);
            auto step_str =  std::to_string(step);
            step_str = std::string(std::to_string(MAX_STEP).length() - step_str.length(), '0') + step_str;
            out.open(step_str+"_waiting_time.txt", std::ofstream::out);
            print_vehicles(out, all_vehicles);
            //img.save(( step_str+prefix_fname).c_str());
            out.close();
        }

        // Stop printing
        /********************************************************************************************************/
        step++;
    }

    MPI_Finalize();
    return 0;
}
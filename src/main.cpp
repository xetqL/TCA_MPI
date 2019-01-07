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

using namespace std;
using namespace tca;

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

unordered_map<long long, Vehicle> parallel_update(const int msx, const int msy,
                                                  const unordered_map<long long,  CA_Cell> &ca_matrix,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map,
                                                  const unordered_map<long long,  Vehicle> &remote_vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new;

    for (const auto &v : vehicles_map) {
        parallel::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, remote_vehicles_map, vehicles_map_new);
    }

    return vehicles_map_new;
}

void generate_vehicles(const int my_rank,
                       const int msx, const int msy, Zoltan_Struct* zz,
                       const unordered_map<long long,  CA_Cell> &ca_matrix,
                             unordered_map<long long,  Vehicle> &vehicles_map) {

    std::vector<std::pair<long long, CA_Cell>> cells;
    std::copy_if(ca_matrix.cbegin(), ca_matrix.cend(), std::back_inserter(cells), [](auto cell){return cell.second.source;});
    //add new vehicles
    std::array<double, 2> pos = {0, 0};
    int PE;
    for (auto &cell : cells) {
        std::tie(pos[0], pos[1]) = cell_to_position(msx, msy, cell.first);
        Zoltan_LB_Point_Assign(zz, &pos.front(), &PE);
        if(my_rank == PE && !exists(vehicles_map, cell.first) && cell.second.has_to_generate(0)) {
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
    if(argc == 5) {
        SIZE_X = std::atoi(argv[3]);
        SIZE_Y = std::atoi(argv[4]);
        ca_matrix = generate_random_manhattan(SIZE_X, SIZE_Y);
    } else {
        std::tie(SIZE_X, SIZE_Y) = read_roadfile(argv[3], &ca_matrix);
    }
    create_random_left_sources(3, SIZE_X, SIZE_Y, &ca_matrix);

    //

    //vector<vector<CA_Cell> > ca_matrix;


    //vector<Vehicle> vehicles;

    unordered_map<long long, Vehicle> vehicle_matrix;
    //vector<vector<Vehicle *> > vehicle_matrix(SIZE_Y, vector<Vehicle *>(SIZE_X, nullptr));

    auto datatype = Vehicle::register_datatype();
    if(!rank)
        randomize_cars_position(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
    auto vehicles = to_vec(vehicle_matrix);
    std::cout << vehicles.size() << std::endl;
    //std::for_each(vehicles.begin(), vehicles.end(), [](auto v){ std::cout << v << std::endl; });

    fprint(out, SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
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

    if(err) std::cout << err << " no esoteric call" << std::endl;

    std::vector<Vehicle> top_vehicles;

    auto zztop = load_balancing::esoteric::divide_data_into_top_bottom(&vehicles , &top_vehicles, incr_cpu, datatype.elements_datatype, bottom);
    //std::for_each(top_vehicles.cbegin(), top_vehicles.cend(), [&rank](auto v){ std::cout << rank << " has " << v << std::endl; });
    auto remote_vehicle = load_balancing::esoteric::exchange(&vehicles, &top_vehicles, zz, zztop, bottom, incr_cpu, datatype.elements_datatype, 1.0);
    //std::for_each(remote_vehicle.cbegin(), remote_vehicle.cend(), [&rank](auto v){ std::cout << rank << " temporary has " << v << std::endl; });
    std::cout << rank << " " << remote_vehicle.size() << std::endl;
    int step = 0;

    while (false) {

        //if(!rank) out.open(prefix_fname + std::to_string(step), std::ofstream::out);
        MPI_Barrier(bottom);

        int recv, sent;
        /*************************************Start parallel exchange********************************************/

        auto remote_data = zoltan_exchange_data(vehicles, zz, datatype.elements_datatype, bottom, recv, sent, 1.2);

        // Stop parallel exchange

        /*************************************Start parallel computation*****************************************/

        auto vehicle_matrix_remote = to_map(SIZE_X, SIZE_Y, remote_data);
        vehicle_matrix = to_map(SIZE_X, SIZE_Y, vehicles);
        vehicle_matrix = parallel_update(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix, vehicle_matrix_remote);
        generate_vehicles(rank, SIZE_X, SIZE_Y, zz, ca_matrix, vehicle_matrix);
        vehicles = to_vec(vehicle_matrix);

        // Stop parallel computation

        /********************************Start load balancing and migration**************************************/

        zoltan_migrate_particles(vehicles, zz, datatype.elements_datatype, bottom);
        zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);

        // Stop load balancing and migration

        /****************************************Start printing**************************************************/

        std::vector<Vehicle> all_vehicles;
        if(wsize > 1)
            gather_elements_on(vehicles, 0, all_vehicles, datatype.elements_datatype, bottom);
        else
            all_vehicles = vehicles;

        if(!rank) {
            auto vehicle_matrix_print = to_map(SIZE_X, SIZE_Y, all_vehicles);
            auto img = zzframe( SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_print);
            auto step_str =  std::to_string(step);
            step_str = std::string(std::to_string(MAX_STEP).length() - step_str.length(), '0') + step_str;
            img.save(( step_str+prefix_fname).c_str());
            //fprint(out, SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_print);
        }

        // Stop printing
        /********************************************************************************************************/
        step++;
        if(!rank) out.close();
    }

    MPI_Finalize();
    return 0;
}
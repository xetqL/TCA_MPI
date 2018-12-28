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

//#include "../include/geometric_load_balancer.hpp"

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
        apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, vehicles_map_new);
    }

    return vehicles_map_new;
}

unordered_map<long long, Vehicle> parallel_update(
                                                  const int msx, const int msy,
                                                  const unordered_map<long long,  CA_Cell> &ca_matrix,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map,
                                                  const unordered_map<long long,  Vehicle> &remote_vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new;
    for (const auto &v : vehicles_map) {
        parallel::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, remote_vehicles_map, vehicles_map_new);
    }
    return vehicles_map_new;
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

    //vector<vector<CA_Cell> > ca_matrix;
    unordered_map<long long, CA_Cell> ca_matrix;

    std::tie(SIZE_X, SIZE_Y) = read_roadfile(argv[1], &ca_matrix);

    //vector<Vehicle> vehicles;

    unordered_map<long long, Vehicle> vehicle_matrix;
    //vector<vector<Vehicle *> > vehicle_matrix(SIZE_Y, vector<Vehicle *>(SIZE_X, nullptr));

    auto datatype = Vehicle::register_datatype();

    randomize_cars_position(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
    auto vehicles = to_vec(vehicle_matrix);

    fprint(out, SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);

    auto zz = zoltan_create_wrapper(ENABLE_AUTOMATIC_MIGRATION, MPI_COMM_WORLD);

    //auto remote_data = zoltan_exchange_data(vehicles, zz, datatype.elements_datatype, bottom, recv, sent);
    //std::for_each(remote_data.begin(), remote_data.end(), [](auto v){std::cout << v << std::endl;});

    if(wsize > 1)  zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);

    int step = 0;
    const std::string prefix_fname(argv[2]);
    while (step < std::atoi(argv[3])) {
        if(!rank) out.open(prefix_fname + std::to_string(step), std::ofstream::out);
        MPI_Barrier(bottom);

        int recv, sent;
        /*************************************Start parallel exchange********************************************/

        auto remote_data = zoltan_exchange_data(vehicles, zz, datatype.elements_datatype, bottom, recv, sent, 1.2);

        // Stop parallel exchange

        /*************************************Start parallel computation*****************************************/

        auto vehicle_matrix_remote = to_map(SIZE_X, SIZE_Y, remote_data);
        vehicle_matrix = to_map(SIZE_X, SIZE_Y, vehicles);
        vehicle_matrix = parallel_update(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix, vehicle_matrix_remote);
        std::cout << vehicle_matrix.size() << std::endl;
        vehicles = to_vec(vehicle_matrix);

        // Stop parallel computation

        /********************************Start load balancing and migration**************************************/

        if(wsize > 1){
            zoltan_migrate_particles(vehicles, zz, datatype.elements_datatype, bottom);
            zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);
        }
        // Stop load balancing and migration

        /****************************************Start printing**************************************************/

        std::vector<Vehicle> all_vehicles;
        if(wsize > 1) gather_elements_on(vehicles, 0, all_vehicles, datatype.elements_datatype, bottom);
        else all_vehicles = vehicles;
        if(!rank) {
            auto vehicle_matrix_print = to_map(SIZE_X, SIZE_Y, all_vehicles);
            fprint(out, SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_print);
        }

        // Stop printing
        /********************************************************************************************************/
        step++;
        if(!rank) out.close();
    }
    MPI_Finalize();
    return 0;
}
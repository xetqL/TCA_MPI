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

int main(int argc, char **argv) {
    int wsize, rank;
    srand(1);
    //MPI_Init(&argc, &argv);
    //MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int SIZE_X, SIZE_Y;

    //vector<vector<CA_Cell> > ca_matrix;
    unordered_map<long long, CA_Cell> ca_matrix;
    std::tie(SIZE_X, SIZE_Y) = read_roadfile("roads.txt", &ca_matrix);

    //vector<Vehicle> vehicles;

    unordered_map<long long, Vehicle> vehicle_matrix;
    //vector<vector<Vehicle *> > vehicle_matrix(SIZE_Y, vector<Vehicle *>(SIZE_X, nullptr));

    //auto com = Vehicle::register_datatype();

    randomize_cars_position(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
    auto vehicles = to_vec(vehicle_matrix);
    std::for_each(vehicles.begin(), vehicles.end(), [](auto v){std::cout << v << std::endl;});
    vehicle_matrix = to_map(SIZE_X, SIZE_Y, vehicles);
    print(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
    //auto zz = zoltan_create_wrapper(ENABLE_AUTOMATIC_MIGRATION, MPI_COMM_WORLD);
    //zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);

    //MPI_Finalize();
    int step = 0;
    const int MAX_NB_VEHICLE = vehicle_matrix.size();

    while (step > -1) {
        usleep(1000000);
        system("clear");
        vehicle_matrix = update(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
        int cars = count_car(SIZE_X, SIZE_Y, vehicle_matrix);
        std::cout << MAX_NB_VEHICLE << " " << count_car(SIZE_X, SIZE_Y, vehicle_matrix) << std::endl;
        print(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
        //assert(MAX_NB_VEHICLE == cars);
        step++;
    }

    return 0;
}
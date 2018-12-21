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

//#include "../include/geometric_load_balancer.hpp"

#define SPATIAL_DISCRETISATION_X = 7.5 //the average length a conventional vehicle occupies in a closely packed jam (and as such, its width is neglected),
#define TEMPORAL_DISCRETISATION = 1.0 //typical driver’s reaction time

enum CA_State {
    Empty, Car, Wall
}; // with int: 0, 1, 2

using namespace std;
using namespace tca;

int char_to_state(char c) {
    switch (c) {
        case 'W':
            return 2;
        case 'C':
            return 1;
        default :
            return 0;
    }
}

std::tuple<std::vector<Vehicle>, vector<vector<Vehicle *> > > update(const int msx, const int msy,
                                                                     const vector<vector<CA_Cell> > &ca_matrix,
                                                                     const std::vector<Vehicle> &vehicles,
                                                                     const vector<vector<Vehicle *> > &vehicles_map) {
    std::vector<Vehicle> vehicles_new;
    vector<vector<Vehicle *> > vehicles_map_new(msy, std::vector<Vehicle *>(msx, nullptr));
    for (const Vehicle &v : vehicles) {
        apply_rule184(msx, msy, ca_matrix, v, vehicles_map, vehicles_new, vehicles_map_new);
    }

    for(auto& v : vehicles_new) {
        int x,y; std::tie(x,y) = v.position;
        std::cout << v << std::endl;
        vehicles_map_new[y][x] = &v;
    }

    return std::make_tuple(vehicles_new, vehicles_map_new);
}

std::tuple<std::vector<Vehicle>, unordered_map<long long, Vehicle *> > update(const int msx, const int msy,
                                                                     const unordered_map<long long,  CA_Cell> &ca_matrix,
                                                                     const std::vector<Vehicle> &vehicles,
                                                                     const unordered_map<long long, Vehicle *> &vehicles_map) {
    std::vector<Vehicle> vehicles_new;
    vehicles_new.reserve(vehicles.size());
    unordered_map<long long, Vehicle *> vehicles_map_new;
    for (const Vehicle &v : vehicles) {
        apply_rule184(msx, msy, ca_matrix, v, vehicles_map, vehicles_new, vehicles_map_new);
    }

    return std::make_tuple(vehicles_new, vehicles_map_new);
}

double speed(int delta_x_meters, int delta_y_seconds) {
    return (delta_x_meters / 1000.0) / (delta_y_seconds / 3600.0);
}

std::pair<int, int> read_roadfile(std::string filename, vector<vector<CA_Cell> > *ca_matrix) {
    std::ifstream f;
    f.open(filename);
    std::string line;
    ca_matrix->clear();
    size_t y = 0;
    while (std::getline(f, line)) {
        size_t x = 0;
        vector<char> v(line.begin(), line.end());
        std::vector<CA_Cell> road;
        for (char c : v) {
            road.emplace_back(std::make_pair(x, y), c);
            x++;
        }
        ca_matrix->push_back(road);
        y++;
    }
    return std::make_pair(ca_matrix->operator[](0).size(), y);
}

std::pair<int, int> read_roadfile(std::string filename, unordered_map<long long, CA_Cell> *ca_matrix) {
    std::ifstream f;
    f.open(filename);
    std::string line;
    ca_matrix->clear();
    size_t y = 0;
    int msx;
    while (std::getline(f, line)) {
        msx = line.size();
        int x = 0;
        vector<char> v(line.begin(), line.end());
        for (char c : v) {
            ca_matrix->operator[](position_to_cell(msx, -1, x, y)) = CA_Cell(c);
            x++;
        }
        y++;
    }
    return std::make_pair(msx, y);
}

void randomize_cars_position(size_t sx, size_t sy, const vector<vector<CA_Cell> > &ca_matrix, vector<Vehicle> *vehicles,
                             vector<vector<Vehicle *>> *vehicle_matrix) {
    int gid = 0;
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix[y][x].direction != NoDirection) {
                if (rand() % 100 > 20) {
                    vehicles->emplace_back(gid, gid, x, y, 1);
                    vehicle_matrix->operator[](y)[x] = &(*(vehicles->end() - 1));
                    gid++;
                }
            }
        }
    }

}

void randomize_cars_position(size_t sx, size_t sy, const unordered_map<long long, CA_Cell> &ca_matrix, vector<Vehicle>& vehicles,
                             unordered_map<long long, Vehicle *>& vehicle_matrix) {

    int gid = 0;
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix.at(position_to_cell(sx, sy, x, y)).direction != NoDirection) {
                if (rand() % 100 > -1) {
                    vehicles.emplace_back(gid, gid, x, y, 1);
                    gid++;
                }
            }
        }
    }

    for(auto& v : vehicles) {
        int x,y; std::tie(x,y) = v.position;
        std::cout << v << std::endl;
        vehicle_matrix[position_to_cell(sx, sy, x, y)] = &v;
    }
}

void randomize_cars_position(int nb_car, size_t sx, size_t sy, const vector<vector<CA_Cell> > &ca_matrix, vector<Vehicle> *vehicles,
                             vector<vector<Vehicle *>> *vehicle_matrix) {
    int gid = 0;
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix[y][x].direction != NoDirection) {
                if (rand() % 100 > 10 && gid < nb_car) {
                    vehicles->emplace_back(gid, gid, x, y, 1);
                    vehicle_matrix->operator[](y)[x] = &(*(vehicles->end() - 1));
                    gid++;
                }
            }
        }
    }
}

void print(size_t sx, size_t sy, const vector<vector<CA_Cell> > &ca_matrix,
           const vector<vector<Vehicle *>> &vehicle_matrix) {
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix.at(y).at(x).direction == NoDirection) std::cout << "*";
            else {
                if (vehicle_matrix.at(y).at(x) == nullptr)
                    std::cout << " ";
                else std::cout << "#";
            }
        }
        std::cout << std::endl;
    }
}

void print(size_t sx, size_t sy, const unordered_map<long long, CA_Cell> &ca_matrix,
           const unordered_map<long long, Vehicle *> &vehicle_matrix) {
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix.at(position_to_cell(sx, sy, x, y)).direction == NoDirection) std::cout << "*";
            else {
                if (get_or_default<long long, Vehicle*>(vehicle_matrix, position_to_cell(sx, sy, x, y), nullptr) == nullptr)
                    std::cout << " ";
                else std::cout << "#";
            }
        }
        std::cout << std::endl;
    }
}

int count_car(int msx, int msy, vector<vector<Vehicle *> >& vehicle_matrix){
    int c = 0;
    for(int i = 0; i < msy; i++){
        for(int j = 0; j < msx; j++){
            if(vehicle_matrix[i][j] != nullptr) c++;
        }
    }
    return c;
}

int main(int argc, char **argv) {
    int wsize, rank;
    srand(time(NULL));
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int SIZE_X, SIZE_Y;

    vector<vector<CA_Cell> > ca_matrix;
    //unordered_map<long long, CA_Cell> ca_matrix;
    std::tie(SIZE_X, SIZE_Y) = read_roadfile("roads.txt", &ca_matrix);

    vector<Vehicle> vehicles;

    //unordered_map<long long, Vehicle *> vehicle_matrix;
    vector<vector<Vehicle *> > vehicle_matrix(SIZE_Y, vector<Vehicle *>(SIZE_X, nullptr));

    auto com = Vehicle::register_datatype();

    randomize_cars_position(SIZE_X, SIZE_Y, ca_matrix, &vehicles, &vehicle_matrix);

    auto zz = zoltan_create_wrapper(ENABLE_AUTOMATIC_MIGRATION, MPI_COMM_WORLD);
    zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);

    MPI_Finalize();

    int step = 0;
    const int MAX_NB_VEHICLE = vehicles.size();
    print(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);

    while (step > -1) {
        usleep(10000);
        system("clear");
        std::tie(vehicles, vehicle_matrix) = update(SIZE_X, SIZE_Y, ca_matrix, vehicles, vehicle_matrix);
        std::cout << MAX_NB_VEHICLE << " " << vehicles.size() << " " << count_car(SIZE_X, SIZE_Y, vehicle_matrix) << std::endl;
        print(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
        assert(MAX_NB_VEHICLE == vehicles.size());
        step++;
    }
    return 0;
}
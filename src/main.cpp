#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <zconf.h>
#include <memory>
#include <map>
#include <cassert>

#define SPATIAL_DISCRETISATION_X = 7.5 //the average length a conventional vehicle occupies in a closely packed jam (and as such, its width is neglected),
#define TEMPORAL_DISCRETISATION = 1.0 //typical driver’s reaction time


enum CA_State {
    Empty, Car, Wall
}; // with int: 0, 1, 2
enum DrivingDirection {
    GoingLeft, GoingRight, GoingUp, GoingDown, Rotary, NoDirection
}; // with int: 0, 1, 2, 3, 4

using namespace std;

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

DrivingDirection char_to_driving(char c) {
    switch (c) {
        case 'L':
            return GoingLeft;
        case 'R':
            return GoingRight;
        case 'U':
            return GoingUp;
        case 'D':
            return GoingDown;
        case 'O':
            return Rotary;
        default :
            return NoDirection;
    }
}

struct CA_Cell {
    DrivingDirection direction;
    std::tuple<int,int> position;
    CA_Cell() : direction(DrivingDirection::NoDirection) {}

    CA_Cell(char c) : direction(char_to_driving(c)) {}
    CA_Cell(std::tuple<int,int> position, char c) : position(std::move(position)), direction(char_to_driving(c)) {}
};

struct Vehicle {
    std::tuple<int, int> position;
    bool rotary_exit_flag;
    float speed;

    Vehicle(int x, int y, float speed) : position(std::make_tuple(x, y)), rotary_exit_flag(rand() % 100 < 30), speed(speed) {}
    Vehicle(std::tuple<int, int> pos, float speed) : position(std::move(pos)), rotary_exit_flag(rand() % 100 < 30), speed(speed) {}
    Vehicle(CA_Cell& c, float speed) : position(c.position), rotary_exit_flag(rand() % 100 < 30), speed(speed) {}
    Vehicle(const CA_Cell& c, float speed) : position(c.position), rotary_exit_flag(rand() % 100 < 30), speed(speed) {}
};

inline void nullptrify(vector<vector<Vehicle *> > &vehicles_map) {
    for (auto &r : vehicles_map)
        for (auto &c : r) c = nullptr;
}
bool can_move(const CA_Cell &destination_cell, const Vehicle* destination_vehicle, const Vehicle* priority_vehicle) {
    return (destination_vehicle == nullptr && destination_cell.direction != Rotary) ||
           (destination_cell.direction == Rotary && (priority_vehicle == nullptr || priority_vehicle->rotary_exit_flag));
}
void in_rotary_rule(const Vehicle* vehicle, const CA_Cell& exit_cell, const CA_Cell& next_rotary_cell, const Vehicle* neighbor_vehicle, std::vector<Vehicle>& vehicles_new, vector<vector<Vehicle *>>& vehicles_map_new){
    int x,y;
    std::tie(x,y) = next_rotary_cell.position;
    if(vehicle->rotary_exit_flag) { // get out, if you can
        if(can_move(exit_cell, neighbor_vehicle, nullptr)) {
            std::tie(x,y) = exit_cell.position;
            vehicles_new.emplace_back(exit_cell, 1);
            vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
        } else {
            vehicles_new.emplace_back(next_rotary_cell, 1);
            vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
        }
    } else {
        vehicles_new.emplace_back(next_rotary_cell, 1);
        vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
    }
}

inline void apply_rule184(const int msx, const int msy,
                   const vector<vector<CA_Cell> > &ca_matrix,
                   const Vehicle &vehicle,
                   const vector<vector<Vehicle *> > &vehicles_map,
                   std::vector<Vehicle>& vehicles_new, vector<vector<Vehicle *>>& vehicles_map_new){

    int x,y;
    std::tie(x, y) = vehicle.position;
    CA_Cell Ni = ca_matrix.at(y).at(x);

    auto D = Ni.direction;
    auto center_cell= ca_matrix.at(y).at(x);
    switch (D) {
        case GoingRight: {
            //if (x+1 == msx) x = -1; //periodic boundary
            if (x + 1 == msx) { // deletion at boundary, if parallel, ask neighbor
                vehicles_map_new[y][x] = nullptr;
            } else {
                auto east = vehicles_map.at(y).at(x + 1);
                auto east_cell = ca_matrix.at(y).at(x + 1);
                if(can_move(east_cell, east, vehicles_map.at(y-1).at(x + 1))) { //
                    vehicles_new.emplace_back(x + 1, y, 1);
                    vehicles_map_new[y][x + 1] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(center_cell, 1);
                    vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
                }
            }
        }
        break;
        case GoingLeft: {
            //if (x+1 == msx) x = -1; //periodic boundary
            if (x == 0) { // deletion at boundary, if parallel, ask neighbor
                vehicles_map_new[y][x] = nullptr;
            } else {
                auto west = vehicles_map.at(y).at(x - 1);
                auto west_cell = ca_matrix.at(y).at(x - 1);
                if(can_move(west_cell, west, vehicles_map.at(y + 1).at(x -1))) { //
                    vehicles_new.emplace_back(x - 1, y, 1);
                    vehicles_map_new[y][x - 1] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(center_cell, 1);
                    vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
                }
            }
        }
            break;
        case GoingUp: {
            //if (x+1 == msx) x = -1; //periodic boundary
            if (y == 0) { // deletion at boundary, if parallel, ask neighbor
                vehicles_map_new[y][x] = nullptr;
            } else {
                auto north = vehicles_map.at(y - 1).at(x);
                auto north_cell = ca_matrix.at(y - 1).at(x);
                if(can_move(north_cell, north, vehicles_map.at(y - 1).at(x - 1))) { //
                    vehicles_new.emplace_back(x, y - 1, 1);
                    vehicles_map_new[y - 1][x] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(center_cell, 1);
                    vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
                }
            }
        }
        break;
        case GoingDown: {
            //if (x+1 == msx) x = -1; //periodic boundary
            if (y + 1 == msy) { // deletion at boundary, if parallel, ask neighbor
                vehicles_map_new[y][x] = nullptr;
            } else {
                auto south = vehicles_map.at(y + 1).at(x);
                auto south_cell = ca_matrix.at(y + 1).at(x);
                if(can_move(south_cell, south, vehicles_map.at(y + 1).at(x + 1))) { //
                    vehicles_new.emplace_back(x, y + 1, 1);
                    vehicles_map_new[y + 1][x] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(center_cell, 1);
                    vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
                }
            }
        }
        break;
        case Rotary: { // von neumann neighborhood
            auto north_cell = ca_matrix.at(y - 1).at(x); auto south_cell = ca_matrix.at(y + 1).at(x);
            auto north_vehicle = vehicles_map.at(y-1).at(x);auto south_vehicle = vehicles_map.at(y+1).at(x);
            auto west_cell  = ca_matrix.at(y).at(x - 1);  auto east_cell  = ca_matrix.at(y).at(x + 1);
            auto west_vehicle = vehicles_map.at(y).at(x - 1);auto east_vehicle = vehicles_map.at(y).at(x + 1);
            //auto vehicle = vehicles_map.at(y).at(x);

            if(south_cell.direction == Rotary && east_cell.direction == Rotary) {   // top left
                in_rotary_rule(&vehicle, west_cell, south_cell, west_vehicle, vehicles_new, vehicles_map_new);
            }
            if(east_cell.direction == Rotary  && north_cell.direction == Rotary){  // bottom left
                in_rotary_rule(&vehicle, south_cell, east_cell, south_vehicle, vehicles_new, vehicles_map_new);
            }
            if(north_cell.direction == Rotary && west_cell.direction == Rotary) {  // bottom right
                in_rotary_rule(&vehicle, east_cell, north_cell, east_vehicle, vehicles_new, vehicles_map_new);
            }
            if(west_cell.direction == Rotary  && south_cell.direction == Rotary) { // top right
                in_rotary_rule(&vehicle, north_cell, west_cell, north_vehicle, vehicles_new, vehicles_map_new);
            }
        }
        break;
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

    return std::make_tuple(vehicles_new, vehicles_map_new);
}

double speed(int delta_x_meters, int delta_y_seconds) {
    return (delta_x_meters / 1000.0) / (delta_y_seconds / 3600.0);
}

std::tuple<int, int> read_roadfile(std::string filename, vector<vector<CA_Cell> > *ca_matrix) {
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
            road.emplace_back(std::make_tuple(x,y), c);
            x++;
        }
        ca_matrix->push_back(road);
        y++;
    }
    return std::make_tuple(ca_matrix->operator[](0).size(), y);
}

void randomize_cars_position(size_t sx, size_t sy, const vector<vector<CA_Cell> > &ca_matrix, vector<Vehicle> *vehicles,
                             vector<vector<Vehicle *>> *vehicle_matrix) {
    int ok = 0;
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix[y][x].direction != NoDirection) {
                if (rand() % 100 > 10) {
                    vehicles->emplace_back(x, y, 1);
                    vehicle_matrix->operator[](y)[x] = &(*(vehicles->end() - 1));
                    ok++;
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

int main(int argc, char **argv) {
    int SIZE_X, SIZE_Y;
    vector<vector<CA_Cell> > ca_matrix;
    std::tie(SIZE_X, SIZE_Y) = read_roadfile("roads.txt", &ca_matrix);

    vector<Vehicle> vehicles;

    vector<vector<Vehicle *> > vehicle_matrix(SIZE_Y, vector<Vehicle *>(SIZE_X, nullptr));

    randomize_cars_position(SIZE_X, SIZE_Y, ca_matrix, &vehicles, &vehicle_matrix);
    print(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);

    std::cout << "Number of vehicles: " << vehicles.size() << std::endl;
    int step = 0;
    while (step >= 0) {
        system("clear");
        int MAX_NB_VEHICLE = vehicles.size();
        std::tie(vehicles, vehicle_matrix) = update(SIZE_X, SIZE_Y, ca_matrix, vehicles, vehicle_matrix);
        assert(MAX_NB_VEHICLE >= vehicles.size());
        print(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
        usleep(200000);
        step++;
    }
    return 0;
}
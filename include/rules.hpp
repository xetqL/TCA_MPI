//
// Created by xetql on 20.12.18.
//

#ifndef CA_ROAD_RULES_HPP
#define CA_ROAD_RULES_HPP


#include <vector>

#include "vehicle.hpp"
#include "ca_cell.hpp"
#include "tca_utils.hpp"

bool can_move(const CA_Cell &destination_cell, const Vehicle* destination_vehicle, const Vehicle* priority_vehicle, bool can_exit_rotary = false) {
    return (destination_cell.direction != Rotary && destination_vehicle == nullptr ) ||
           (destination_cell.direction == Rotary &&
                ((priority_vehicle  == nullptr && destination_vehicle == nullptr) ||
                 (priority_vehicle  != nullptr && (priority_vehicle->rotary_exit_flag && can_exit_rotary) && destination_vehicle == nullptr)) );
}

std::pair<std::pair<int, int>, std::pair<int, int>> get_rotary_neighbors(int x, int y,
                                                                         const std::vector<std::vector<CA_Cell> > &ca_matrix) {
    auto north_cell = ca_matrix.at(y - 1).at(x);  auto south_cell = ca_matrix.at(y + 1).at(x);
    auto west_cell  = ca_matrix.at(y).at(x - 1);  auto east_cell  = ca_matrix.at(y).at(x + 1);

    if(south_cell.direction == Rotary && east_cell.direction == Rotary) {   // top left
        return std::make_pair(west_cell.position, south_cell.position);
    }
    if(east_cell.direction == Rotary  && north_cell.direction == Rotary){  // bottom left
        return std::make_pair(south_cell.position, east_cell.position);
    }
    if(north_cell.direction == Rotary && west_cell.direction == Rotary) {  // bottom right
        return std::make_pair(east_cell.position, north_cell.position);
    }
    if(west_cell.direction == Rotary  && south_cell.direction == Rotary) { // top right
        return std::make_pair(north_cell.position, west_cell.position);
    }
}

std::pair<std::pair<int, int>, std::pair<int, int>> get_rotary_neighbors(int msx, int msy, int x, int y,
                                                                         const std::unordered_map<long long, CA_Cell> &ca_matrix) {
    long long xy;

    xy = position_to_cell(msx, msy, x, y - 1 == -1 ? msy : y - 1); auto north_cell = get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());
    xy = position_to_cell(msx, msy, x, y + 1 == msy ? -1 : y + 1); auto south_cell = get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());
    xy = position_to_cell(msx, msy, x - 1 == -1 ? msx : x - 1, y); auto west_cell  = get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());
    xy = position_to_cell(msx, msy, x + 1 == msx ? -1 : x + 1, y); auto east_cell  = get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());

    if(south_cell.direction == Rotary && east_cell.direction == Rotary) {   // top left
        return std::make_pair(west_cell.position, south_cell.position);
    }
    if(east_cell.direction == Rotary  && north_cell.direction == Rotary){  // bottom left
        return std::make_pair(south_cell.position, east_cell.position);
    }
    if(north_cell.direction == Rotary && west_cell.direction == Rotary) {  // bottom right
        return std::make_pair(east_cell.position, north_cell.position);
    }
    if(west_cell.direction == Rotary  && south_cell.direction == Rotary) { // top right
        return std::make_pair(north_cell.position, west_cell.position);
    }
}
namespace deprecated {
void in_rotary_rule(const Vehicle* vehicle,
                    const CA_Cell& exit_cell,
                    const CA_Cell& next_rotary_cell,
                    const Vehicle* neighbor_vehicle,
                    std::vector<Vehicle>& vehicles_new,
                    std::vector<std::vector<Vehicle *>>& vehicles_map_new){
    int x,y;
    std::tie(x,y) = next_rotary_cell.position;
    if(vehicle->rotary_exit_flag) { // get out, if you can
        if(can_move(exit_cell, neighbor_vehicle, nullptr)) {
            std::tie(x,y) = exit_cell.position;
            vehicles_new.emplace_back(vehicle->gid, vehicle->lid, exit_cell, 1);
            //vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
        } else {
            vehicles_new.emplace_back(vehicle->gid, vehicle->lid, next_rotary_cell, 1);
            //vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
        }
    } else {
        vehicles_new.emplace_back(vehicle->gid, vehicle->lid, next_rotary_cell, 1);
        //vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
    }
}

void in_road_rule(const Vehicle* vehicle,
                  const CA_Cell& next_cell,
                  const Vehicle* neighbor_vehicle,
                  const Vehicle* priority_vehicle,
                  const std::vector<std::vector<CA_Cell> > &ca_matrix,
                  const std::vector<std::vector<Vehicle *> > &vehicles_map,
                  std::vector<Vehicle>& vehicles_new,
                  std::vector<std::vector<Vehicle *>>& vehicles_map_new) {
    int x, y;
    bool priority_vehicle_can_exit_rotary = false;
    if(next_cell.direction == Rotary && priority_vehicle != nullptr){
        std::tie(x,y) = priority_vehicle->position;
        std::pair<int,int> exit, next; std::tie(exit, next) = get_rotary_neighbors(x, y, ca_matrix);
        //priority_vehicle_can_exit_rotary = vehicles_map[exit.second][exit.first] == nullptr;
    }

    if(can_move(next_cell, neighbor_vehicle, priority_vehicle, priority_vehicle_can_exit_rotary)) { //
        std::tie(x, y) = next_cell.position;
    }else {
        std::tie(x, y) = vehicle->position;
    }
    vehicles_new.emplace_back(vehicle->gid, vehicle->lid, x, y, 1);
    //vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
}


inline void apply_rule184(const int msx, const int msy,
                          const std::vector<std::vector<CA_Cell> > &ca_matrix,
                          const Vehicle &vehicle,
                          const std::vector<std::vector<Vehicle *> > &vehicles_map,
                          std::vector<Vehicle>& vehicles_new, std::vector<std::vector<Vehicle *>>& vehicles_map_new){
    int x, y;
    std::tie(x, y) = vehicle.position;
    CA_Cell Ni = ca_matrix.at(y).at(x);
    auto D = Ni.direction;
    switch (D) {
        case GoingRight: {
            if (x + 1 == msx) { // deletion at boundary, if parallel, ask neighbor
                // vehicles_map_new[y][x] = nullptr;
                x = -1;
            }   // else
            {
                auto next_cell = ca_matrix.at(y).at(x + 1); auto next_vehicle = vehicles_map.at(y).at(x + 1);
                auto priority_vehicle = vehicles_map.at(y - 1).at(x + 1);
                in_road_rule(&vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_new, vehicles_map_new);
            }
        }
            break;
        case GoingLeft: {
            if (x == 0) { // deletion at boundary, if parallel, ask neighbor
                // vehicles_map_new[y][x] = nullptr;
                x = msx;
            }   // else
            {
                auto next_vehicle = vehicles_map.at(y).at(x - 1); auto next_cell = ca_matrix.at(y).at(x - 1);
                auto priority_vehicle = vehicles_map.at(y + 1).at(x -1);
                in_road_rule(&vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map,vehicles_new, vehicles_map_new);
            }
        }
            break;
        case GoingUp: {
            if (y == 0) { // deletion at boundary, if parallel, ask neighbor
                // vehicles_map_new[y][x] = nullptr;
                y = msy;
            }   // else
            {
                auto next_vehicle = vehicles_map.at(y - 1).at(x); auto next_cell = ca_matrix.at(y - 1).at(x);
                auto priority_vehicle = vehicles_map.at(y - 1).at(x - 1);
                in_road_rule(&vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_new, vehicles_map_new);
            }
        }
            break;
        case GoingDown: {
            if (y + 1 == msy) { // deletion at boundary, if parallel, ask neighbor
                // vehicles_map_new[y][x] = nullptr;
                y = -1;
            }   // else
            {
                auto next_vehicle = vehicles_map.at(y + 1).at(x); auto next_cell = ca_matrix.at(y + 1).at(x);
                auto priority_vehicle = vehicles_map.at(y + 1).at(x + 1);
                in_road_rule(&vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_new, vehicles_map_new);
            }
        }
            break;
        case Rotary: { // von neumann neighborhood
            std::pair<int,int> exit, next; std::tie(exit, next) = get_rotary_neighbors(x, y, ca_matrix);

            auto exit_cell    = ca_matrix.at(exit.second).at(exit.first);
            auto next_cell    = ca_matrix.at(next.second).at(next.first);
            auto exit_vehicle = vehicles_map.at(exit.second).at(exit.first);

            in_rotary_rule(&vehicle, exit_cell, next_cell, exit_vehicle, vehicles_new, vehicles_map_new);

        }
            break;
    }
}

}
namespace sequential {
void in_rotary_rule(const int msx, const int msy,
                    const Vehicle* vehicle,
                    const CA_Cell& exit_cell,
                    const CA_Cell& next_rotary_cell,
                    const Vehicle* neighbor_vehicle,
                    std::unordered_map<long long, Vehicle>& vehicles_map_new){
    int x,y;
    std::tie(x,y) = next_rotary_cell.position;
    long long xy = position_to_cell(msx, msy, x, y);
    //std::cout << vehicle->rotary_exit_flag;
    if(vehicle->rotary_exit_flag) { // get out, if you can
        if(can_move(exit_cell, neighbor_vehicle, nullptr)) {
            std::tie(x,y) = exit_cell.position;
            xy = position_to_cell(msx, msy, x, y);
            vehicles_map_new[xy] = Vehicle(vehicle->gid, vehicle->lid, exit_cell, 1);
            return;
        }
    }
    vehicles_map_new[xy] = Vehicle(vehicle->gid, vehicle->lid, next_rotary_cell, 1);
}

void in_road_rule(const int msx, const int msy,
                  const Vehicle* vehicle,
                  const CA_Cell& next_cell,
                  const Vehicle* neighbor_vehicle,
                  const Vehicle* priority_vehicle,
                  const std::unordered_map<long long, CA_Cell> &ca_matrix,
                  const std::unordered_map<long long, Vehicle> &vehicles_map,
                  std::unordered_map<long long, Vehicle> &vehicles_map_new) {
    int x, y;
    bool priority_vehicle_can_exit_rotary = false;

    if(next_cell.direction == Rotary && priority_vehicle != nullptr){
        std::tie(x,y) = priority_vehicle->position;
        std::pair<int,int> exit, next; std::tie(exit, next) = get_rotary_neighbors(msx, msy, x, y, ca_matrix);
        priority_vehicle_can_exit_rotary = !exists(vehicles_map, position_to_cell(msx, msy, exit));
        //priority_vehicle_can_exit_rotary = !exists(vehicles_map, position_to_cell(msx, msy, next));

    }
    std::pair<int,int> p;
    if(can_move(next_cell, neighbor_vehicle, priority_vehicle, priority_vehicle_can_exit_rotary)) { //
        p = next_cell.position;
    }else {
        p = vehicle->position;
    }
    std::tie(x,y) = p;
    // << "From (" << vehicle->position.first << "," << vehicle->position.second << ") to " << "(" << x << "," << y << ")" << std::endl;
    long long xy = position_to_cell(msx, msy, x, y);
    vehicles_map_new[xy] = Vehicle(vehicle->gid, vehicle->lid, x, y, 1);
}

inline void apply_rule184(const int msx, const int msy,
                          const std::unordered_map<long long, CA_Cell> &ca_matrix,
                          const Vehicle &vehicle,
                          const std::unordered_map<long long, Vehicle> &vehicles_map,
                          std::unordered_map<long long, Vehicle> &vehicles_map_new) {
    int x, y;
    std::tie(x, y) = vehicle.position;
    long long xy = position_to_cell(msx, msy, x, y);
    CA_Cell Ni = ca_matrix.at(xy);
    auto D = Ni.direction;

    switch (D) {
        case GoingRight: {
            if (x + 1 == msx) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                break;
                //x = -1;
            }
            xy = position_to_cell(msx, msy, x+1, y);
            auto next_cell = ca_matrix.at(xy); //will always be ok.
            auto next_vehicle = exists(vehicles_map, xy) ? &vehicles_map.at(xy) : nullptr;

            xy = position_to_cell(msx, msy, x+1, y-1);
            auto priority_vehicle = exists(vehicles_map, xy) && ca_matrix.at(xy).direction == Rotary  ? &vehicles_map.at(xy) : nullptr;

            in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix,  vehicles_map, vehicles_map_new);

        }
            break;
        case GoingLeft: {

            if (x == 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                x = msx;
            }

            xy = position_to_cell(msx, msy, x-1, y);
            auto next_cell = ca_matrix.at(xy);
            auto next_vehicle = exists(vehicles_map, xy) ? &vehicles_map.at(xy) : nullptr;
            xy = position_to_cell(msx, msy, x-1, y+1);
            auto priority_vehicle = exists(vehicles_map, xy) && ca_matrix.at(xy).direction == Rotary ? &vehicles_map.at(xy) : nullptr;
            in_road_rule(msx,msy,&vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map,  vehicles_map_new);

        }
            break;
        case GoingUp: {

            if (y == 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;

                y = msy;
            }

            xy = position_to_cell(msx, msy, x, y - 1 == -1 ? msy : y - 1);
            auto next_cell = ca_matrix.at(xy);
            auto next_vehicle = exists(vehicles_map, xy) ? &vehicles_map.at(xy) : nullptr;
            xy = position_to_cell(msx, msy, x-1 == -1 ? msx : x - 1, y - 1 == -1 ? msy : y - 1);
            auto priority_vehicle = exists(vehicles_map, xy) && ca_matrix.at(xy).direction == Rotary  ? &vehicles_map.at(xy) : nullptr;
            in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix,  vehicles_map, vehicles_map_new);

        }
            break;
        case GoingDown: {
            if(y+1 == msy) break;

            xy = position_to_cell(msx, msy, x, y + 1 == msy ? 0 : y+1);
            auto next_cell = ca_matrix.at(xy);
            auto next_vehicle = exists(vehicles_map, xy) ? &vehicles_map.at(xy) : nullptr;
            xy = position_to_cell(msx, msy, x+1 == msx ? 0 : x+1, y+1 == msy ? 0 : y+1);
            auto priority_vehicle = exists(vehicles_map, xy) && ca_matrix.at(xy).direction == Rotary ? &vehicles_map.at(xy) : nullptr;
            in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix,  vehicles_map,vehicles_map_new);

        }
            break;
        case Rotary: { // von neumann neighborhood

            std::pair<int,int> exit, next; std::tie(exit, next) = get_rotary_neighbors(msx, msy, x, y, ca_matrix);
            auto exit_cell    = ca_matrix.at(position_to_cell(msx, msy, exit));
            auto next_cell    = ca_matrix.at(position_to_cell(msx, msy, next));
            xy = position_to_cell(msx,msy, exit);
            auto exit_vehicle = exists(vehicles_map, xy) ? &vehicles_map.at(xy) : nullptr;

            in_rotary_rule(msx, msy, &vehicle, exit_cell, next_cell, exit_vehicle, vehicles_map_new);
            //// << " Rotary" << position_to_cell(msx,msy,x,y) << std::endl;

        }
            break;
    }
}
}
namespace parallel {
void in_road_rule(const int msx, const int msy,
                  const Vehicle* vehicle,
                  const CA_Cell& next_cell,
                  const Vehicle* neighbor_vehicle,
                  const Vehicle* priority_vehicle,
                  const std::unordered_map<long long, CA_Cell> &ca_matrix,
                  const std::unordered_map<long long, Vehicle> &vehicles_map,
                  const std::unordered_map<long long, Vehicle> &vehicles_map_remote,
                  std::unordered_map<long long, Vehicle> &vehicles_map_new){
    int x, y;
    bool priority_vehicle_can_exit_rotary = false;

    if(next_cell.direction == Rotary && priority_vehicle != nullptr){
        std::tie(x,y) = priority_vehicle->position;
        std::pair<int,int> exit, next; std::tie(exit, next) = get_rotary_neighbors(msx, msy, x, y, ca_matrix);
        priority_vehicle_can_exit_rotary =
                !exists(vehicles_map, position_to_cell(msx, msy, exit)) &&
                !exists(vehicles_map_remote, position_to_cell(msx, msy, exit));
        //priority_vehicle_can_exit_rotary = !exists(vehicles_map, position_to_cell(msx, msy, next));

    }
    std::pair<int,int> p;
    if(can_move(next_cell, neighbor_vehicle, priority_vehicle, priority_vehicle_can_exit_rotary)) { //
        p = next_cell.position;
    }else {
        p = vehicle->position;
    }
    std::tie(x,y) = p;
    // << "From (" << vehicle->position.first << "," << vehicle->position.second << ") to " << "(" << x << "," << y << ")" << std::endl;
    long long xy = position_to_cell(msx, msy, x, y);
    vehicles_map_new[xy] = Vehicle(vehicle->gid, vehicle->lid, x, y, 1);
}

void apply_rule184(const int msx, const int msy,
                          const std::unordered_map<long long, CA_Cell> &ca_matrix,
                          const Vehicle &vehicle,
                          const std::unordered_map<long long, Vehicle> &vehicles_map,
                          const std::unordered_map<long long, Vehicle> &vehicles_map_remote,
                          std::unordered_map<long long, Vehicle> &vehicles_map_new){
    int x, y;
    std::tie(x, y) = vehicle.position;
    long long xy = position_to_cell(msx, msy, x, y);
    CA_Cell Ni = ca_matrix.at(xy);
    auto D = Ni.direction;
    switch (D) {
        case GoingRight: {

            if (x + 1 == msx) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                break;
                //x = -1;
            }
            xy = position_to_cell(msx, msy, x+1, y);
            auto next_cell = ca_matrix.at(xy); //will always be ok.
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            //auto next_vehicle = next_vehicle_iterator != vehicles_map_remote.end()? : &(*next_vehicle_iterator) : nullptr;
            xy = position_to_cell(msx, msy, x+1, y-1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](auto v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;

            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case GoingLeft: {

            if (x == 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;

                break; x = msx;
            }

            xy = position_to_cell(msx, msy, x-1, y);
            auto next_cell = ca_matrix.at(xy);
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            xy = position_to_cell(msx, msy, x-1, y+1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](auto v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;
            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case GoingUp: {

            if (y == 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                break; y = msy;
            }

            xy = position_to_cell(msx, msy, x, y - 1 == -1 ? msy : y - 1);
            auto next_cell = ca_matrix.at(xy);
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            xy = position_to_cell(msx, msy, x-1 == -1 ? msx : x - 1, y - 1 == -1 ? msy : y - 1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](auto v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;
            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case GoingDown: {
            if(y+1 == msy) break;

            xy = position_to_cell(msx, msy, x, y + 1 == msy ? 0 : y+1);
            auto next_cell = ca_matrix.at(xy);
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            xy = position_to_cell(msx, msy, x+1 == msx ? 0 : x+1, y+1 == msy ? 0 : y+1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](auto v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;
            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case Rotary: { // von neumann neighborhood

            std::pair<int,int> exit, next; std::tie(exit, next) = get_rotary_neighbors(msx, msy, x, y, ca_matrix);
            auto exit_cell    = ca_matrix.at(position_to_cell(msx, msy, exit));
            auto next_cell    = ca_matrix.at(position_to_cell(msx, msy, next));
            xy = position_to_cell(msx,msy, exit);
            const Vehicle* exit_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);

            sequential::in_rotary_rule(msx, msy, &vehicle, exit_cell, next_cell, exit_vehicle, vehicles_map_new);
            //// << " Rotary" << position_to_cell(msx,msy,x,y) << std::endl;

        }
            break;
        case NoDirection:
            std::cerr << "No direction: skipping..." << std::endl;
    }
}
}


#endif //CA_ROAD_RULES_HPP

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
           (priority_vehicle  == nullptr || (priority_vehicle->rotary_exit_flag && can_exit_rotary)));
}
std::pair<std::pair<int, int>, std::pair<int, int>> get_rotary_neighbors(int x, int y,
                                                                         const std::vector<std::vector<CA_Cell> > &ca_matrix) {
    auto north_cell = ca_matrix.at(y - 1).at(x); auto south_cell = ca_matrix.at(y + 1).at(x);
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

void in_rotary_rule(const int msx, const int msy,
                    const Vehicle* vehicle,
                    const CA_Cell& exit_cell,
                    const CA_Cell& next_rotary_cell,
                    const Vehicle* neighbor_vehicle,
                    std::vector<Vehicle>& vehicles_new,
                    std::unordered_map<long long, Vehicle *>& vehicles_map_new){
    int x,y;
    std::tie(x,y) = next_rotary_cell.position;
    long long xy = position_to_cell(msx, msy, x, y);
    if(vehicle->rotary_exit_flag) { // get out, if you can
        if(can_move(exit_cell, neighbor_vehicle, nullptr)) {
            std::tie(x,y) = exit_cell.position;
            xy = position_to_cell(msx, msy, x, y);
            vehicles_new.emplace_back(vehicle->gid, vehicle->lid, exit_cell, 1);
            //vehicles_map_new[xy] = &(*(vehicles_new.end() - 1));
        } else {
            vehicles_new.emplace_back(vehicle->gid, vehicle->lid, next_rotary_cell, 1);
            //vehicles_map_new[xy] = &(*(vehicles_new.end() - 1));
        }
    } else {
        vehicles_new.emplace_back(vehicle->gid, vehicle->lid, next_rotary_cell, 1);
        //vehicles_map_new[xy] = &(*(vehicles_new.end() - 1));
    }
}

void in_road_rule(const int msx, const int msy,
                  const Vehicle* vehicle,
                  const CA_Cell& next_cell,
                  const Vehicle* neighbor_vehicle,
                  const Vehicle* priority_vehicle,
                  std::vector<Vehicle>& vehicles_new,
                  std::unordered_map<long long, Vehicle *>& vehicles_map_new) {
    int x, y;

    if(can_move(next_cell, neighbor_vehicle, priority_vehicle)) { //
        std::tie(x, y) = next_cell.position;
    }else {
        std::tie(x, y) = vehicle->position;
    }
    long long xy = position_to_cell(msx, msy, x, y);
    vehicles_new.emplace_back(vehicle->gid, vehicle->lid, x, y, 1);
    //vehicles_map_new[xy] = &(*(vehicles_new.end() - 1));
    //std::cout << vehicles_map_new[y].size() << std::endl;
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
        std::cout << x << std::endl;
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

            auto exit_cell = ca_matrix.at(exit.second).at(exit.first);
            auto next_cell = ca_matrix.at(next.second).at(next.first);
            auto exit_vehicle = vehicles_map.at(exit.second).at(exit.first);

            in_rotary_rule(&vehicle, exit_cell, next_cell, exit_vehicle, vehicles_new, vehicles_map_new);

        }
            break;
    }
}

inline void apply_rule184(const int msx, const int msy,
                          const std::unordered_map<long long, CA_Cell> &ca_matrix,
                          const Vehicle &vehicle,
                          const std::unordered_map<long long, Vehicle *> &vehicles_map,
                          std::vector<Vehicle>& vehicles_new,
                          std::unordered_map<long long, Vehicle *> &vehicles_map_new) {
    int x, y;
    std::tie(x, y) = vehicle.position;
    long long xy = position_to_cell(msx, msy, x, y);
    CA_Cell Ni = ca_matrix.at(xy);
    auto D = Ni.direction;
    switch (D) {
        case GoingRight: {
            if (x + 1 == msx) { // deletion at boundary, if parallel, ask neighbor
            //    vehicles_map_new[y][x] = nullptr;
                x = -1;
            }
            {
                xy = position_to_cell(msx, msy, x+1, y);
                auto next_cell = ca_matrix.at(xy); //will always be ok.
                auto next_vehicle = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
                xy = position_to_cell(msx, msy, x+1, y-1);

                auto priority_vehicle = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
                in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, vehicles_new, vehicles_map_new);
            }
        }
            break;
        case GoingLeft: {
            if (x == 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                x = msx;
            }
            {
                xy = position_to_cell(msx, msy, x-1, y);
                auto next_cell = ca_matrix.at(xy);
                auto next_vehicle = get_or_default<long long, Vehicle*>(vehicles_map,xy, nullptr);
                xy = position_to_cell(msx, msy, x-1, y+1);
                auto priority_vehicle = get_or_default<long long , Vehicle*>(vehicles_map,xy, nullptr);
                in_road_rule(msx,msy,&vehicle, next_cell, next_vehicle, priority_vehicle, vehicles_new, vehicles_map_new);
            }
        }
            break;
        case GoingUp: {
            if (y == 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                y = msy;
            }
            {
                xy = position_to_cell(msx, msy, x, y-1);
                auto next_cell = ca_matrix.at(xy);
                auto next_vehicle = get_or_default<long long, Vehicle*>(vehicles_map,xy, nullptr);
                xy = position_to_cell(msx, msy, x-1, y-1);
                auto priority_vehicle = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
                in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, vehicles_new, vehicles_map_new);
            }
        }
            break;
        case GoingDown: {
            if (y + 1 == msy) { // deletion at boundary, if parallel, ask neighbor //vehicles_map_new[y][x] = nullptr;
                y = -1;
            }
            {
                xy = position_to_cell(msx, msy, x, y+1);
                auto next_cell = ca_matrix.at(xy);
                auto next_vehicle = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
                xy = position_to_cell(msx, msy, x+1, y+1);
                auto priority_vehicle = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
                in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, vehicles_new, vehicles_map_new);
            }
        }
            break;
        case Rotary: { // von neumann neighborhood
            xy = position_to_cell(msx, msy, x, y-1); auto north_cell = ca_matrix.at(xy); auto north_vehicle = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
            xy = position_to_cell(msx, msy, x, y+1); auto south_cell = ca_matrix.at(xy); auto south_vehicle = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
            xy = position_to_cell(msx, msy, x-1, y); auto west_cell  = ca_matrix.at(xy); auto west_vehicle  = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);
            xy = position_to_cell(msx, msy, x+1, y); auto east_cell  = ca_matrix.at(xy); auto east_vehicle  = get_or_default<long long, Vehicle*>(vehicles_map, xy, nullptr);

            if(south_cell.direction == Rotary && east_cell.direction == Rotary) {   // top left
                in_rotary_rule(msx, msy, &vehicle, west_cell, south_cell, west_vehicle, vehicles_new, vehicles_map_new);
            }
            if(east_cell.direction == Rotary  && north_cell.direction == Rotary){  // bottom left
                in_rotary_rule(msx, msy, &vehicle, south_cell, east_cell, south_vehicle, vehicles_new, vehicles_map_new);
            }
            if(north_cell.direction == Rotary && west_cell.direction == Rotary) {  // bottom right
                in_rotary_rule(msx, msy, &vehicle, east_cell, north_cell, east_vehicle, vehicles_new, vehicles_map_new);
            }
            if(west_cell.direction == Rotary  && south_cell.direction == Rotary) { // top right
                in_rotary_rule(msx, msy, &vehicle, north_cell, west_cell, north_vehicle, vehicles_new, vehicles_map_new);
            }
        }
            break;
    }
}

#endif //CA_ROAD_RULES_HPP

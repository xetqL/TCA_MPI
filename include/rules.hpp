//
// Created by xetql on 20.12.18.
//

#ifndef CA_ROAD_RULES_HPP
#define CA_ROAD_RULES_HPP


#include <vector>

#include "vehicle.hpp"
#include "ca_cell.hpp"

bool can_move(const CA_Cell &destination_cell, const Vehicle* destination_vehicle, Vehicle* priority_vehicle) {
    return (destination_vehicle == nullptr && destination_cell.direction != Rotary) ||
           (destination_cell.direction == Rotary && (priority_vehicle == nullptr || priority_vehicle->rotary_exit_flag));
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
            vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
        } else {
            vehicles_new.emplace_back(vehicle->gid, vehicle->lid, next_rotary_cell, 1);
            vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
        }
    } else {
        vehicles_new.emplace_back(vehicle->gid, vehicle->lid, next_rotary_cell, 1);
        vehicles_map_new[y][x] = &(*(vehicles_new.end() - 1));
    }
}

inline void apply_rule184(const int msx, const int msy,
                          const std::vector<std::vector<CA_Cell> > &ca_matrix,
                          const Vehicle &vehicle,
                          const std::vector<std::vector<Vehicle *> > &vehicles_map,
                          std::vector<Vehicle>& vehicles_new, std::vector<std::vector<Vehicle *>>& vehicles_map_new){

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
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, x + 1, y, 1);
                    vehicles_map_new[y][x + 1] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, center_cell, 1);
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
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, x - 1, y, 1);
                    vehicles_map_new[y][x - 1] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, center_cell, 1);
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
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, x, y - 1, 1);
                    vehicles_map_new[y - 1][x] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, center_cell, 1);
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
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, x, y + 1, 1);
                    vehicles_map_new[y + 1][x] = &(*(vehicles_new.end() - 1));
                }else {
                    vehicles_new.emplace_back(vehicle.gid, vehicle.lid, center_cell, 1);
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

#endif //CA_ROAD_RULES_HPP

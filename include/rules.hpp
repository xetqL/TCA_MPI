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

std::tuple<int, std::pair<int, int>, std::pair<int, int>> get_rotary_neighbors(int msx, int msy, int x, int y,
                                                                         const std::unordered_map<long long, CA_Cell>& ca_matrix) {
    long long xy;
    const long long max_xy = msx*msy;
    xy = position_to_cell(msx, msy, x, y - 1 == -1 ? msy : y - 1); auto north_cell = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix.at(xy);
    xy = position_to_cell(msx, msy, x, y + 1 == msy ? -1 : y + 1); auto south_cell = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix.at(xy); //get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());
    xy = position_to_cell(msx, msy, x - 1 == -1 ? msx : x - 1, y); auto west_cell  = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix.at(xy);//get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());
    xy = position_to_cell(msx, msy, x + 1 == msx ? -1 : x + 1, y); auto east_cell  = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix.at(xy);//get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());

    if(south_cell.direction == Rotary && east_cell.direction == Rotary) {   // top left
        return std::make_tuple(1, west_cell.position, south_cell.position);
    }
    if(east_cell.direction == Rotary  && north_cell.direction == Rotary){  // bottom left
        return std::make_tuple(2, south_cell.position, east_cell.position);
    }
    if(north_cell.direction == Rotary && west_cell.direction == Rotary) {  // bottom right
        return std::make_tuple(3, east_cell.position, north_cell.position);
    }
    if(west_cell.direction == Rotary  && south_cell.direction == Rotary) { // top right
        return std::make_tuple(4, north_cell.position, west_cell.position);
    }
}
std::tuple<int, std::pair<int, int>, std::pair<int, int>> get_rotary_neighbors(int msx, int msy, int x, int y,
                                                                         const std::unordered_map<long long, CA_Cell>* ca_matrix) {
    long long xy;
    const long long max_xy = msx*msy;
    xy = position_to_cell(msx, msy, x, y - 1 == -1 ? msy : y - 1); auto north_cell = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix->at(xy);
    xy = position_to_cell(msx, msy, x, y + 1 == msy ? -1 : y + 1); auto south_cell = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix->at(xy); //get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());
    xy = position_to_cell(msx, msy, x - 1 == -1 ? msx : x - 1, y); auto west_cell  = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix->at(xy);//get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());
    xy = position_to_cell(msx, msy, x + 1 == msx ? -1 : x + 1, y); auto east_cell  = xy >= max_xy || xy < 0? CA_Cell::get_fake_cell() : ca_matrix->at(xy);//get_or_default(ca_matrix, xy, CA_Cell::get_fake_cell());

    if(south_cell.direction == Rotary && east_cell.direction == Rotary) {   // top left
        return std::make_tuple(1, west_cell.position, south_cell.position);
    }
    if(east_cell.direction == Rotary  && north_cell.direction == Rotary){  // bottom left
        return std::make_tuple(2, south_cell.position, east_cell.position);
    }
    if(north_cell.direction == Rotary && west_cell.direction == Rotary) {  // bottom right
        return std::make_tuple(3, east_cell.position, north_cell.position);
    }
    if(west_cell.direction == Rotary  && south_cell.direction == Rotary) { // top right
        return std::make_tuple(4, north_cell.position, west_cell.position);
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
            vehicles_map_new[xy] = vehicle->drive(exit_cell.position, 1);
            return;
        }
    }
    vehicles_map_new[xy] = vehicle->drive(next_rotary_cell.position, 1);
    vehicles_map_new[xy].waiting_time++;
}

void in_road_rule(const int msx, const int msy,
                  const Vehicle* vehicle,
                  const CA_Cell& next_cell,
                  const Vehicle* neighbor_vehicle,
                  const Vehicle* priority_vehicle,
                  const std::unordered_map<long long, CA_Cell> &ca_matrix,
                  const std::unordered_map<long long, Vehicle> &vehicles_map,
                  std::unordered_map<long long, Vehicle> &vehicles_map_new) {
    int x, y, rid;
    bool priority_vehicle_can_exit_rotary = false;

    if(next_cell.direction == Rotary && priority_vehicle != nullptr){
        std::tie(x,y) = priority_vehicle->position;
        std::pair<int,int> exit, next; std::tie(rid, exit, next) = get_rotary_neighbors(msx, msy, x, y, ca_matrix);
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
    int x, y, rid;
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

            std::pair<int,int> exit, next; std::tie(rid, exit, next) = get_rotary_neighbors(msx, msy, x, y, ca_matrix);
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
} //end of namespace: sequential
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
    int x, y, rid;
    bool priority_vehicle_can_exit_rotary = false;

    if(next_cell.direction == Rotary && priority_vehicle != nullptr){
        std::tie(x,y) = priority_vehicle->position;
        std::pair<int,int> exit, next; std::tie(rid, exit, next) = get_rotary_neighbors(msx, msy, x, y, &ca_matrix);
        priority_vehicle_can_exit_rotary =
                !exists(vehicles_map, position_to_cell(msx, msy, exit)) &&
                !exists(vehicles_map_remote, position_to_cell(msx, msy, exit));
        //priority_vehicle_can_exit_rotary = !exists(vehicles_map, position_to_cell(msx, msy, next));
    }

    std::pair<int,int> p;
    int waiting_time = 0;
    std::tie(x, y) = vehicle->position;
    long long xy = position_to_cell(msx, msy, x, y);
    float speed = vehicle->speed;
    if(speed <= 0 || !can_move(next_cell, neighbor_vehicle, priority_vehicle, priority_vehicle_can_exit_rotary)){
        p = vehicle->position;
        waiting_time = vehicle->waiting_time+1;
    } else p = next_cell.position;

    std::tie(x, y) = p;
    // << "From (" << vehicle->position.first << "," << vehicle->position.second << ") to " << "(" << x << "," << y << ")" << std::endl;
    xy = position_to_cell(msx, msy, x, y);
    vehicles_map_new[xy] = vehicle->drive(x, y, speed);
    vehicles_map_new[xy].waiting_time = waiting_time;
}

void apply_rule184(const int msx, const int msy,
                          const std::unordered_map<long long, CA_Cell> &ca_matrix,
                          const Vehicle &vehicle,
                          const std::unordered_map<long long, Vehicle> &vehicles_map,
                          const std::unordered_map<long long, Vehicle> &vehicles_map_remote,
                          std::unordered_map<long long, Vehicle> *_vehicles_map_new){
    int x, y, rid;
    std::unordered_map<long long, Vehicle>& vehicles_map_new = *_vehicles_map_new;
    std::tie(x, y) = vehicle.position;
    long long xy = position_to_cell(msx, msy, x, y);
    CA_Cell Ni = ca_matrix.at(xy);
    auto D = Ni.direction;

    //std::cout << "Speed: "<< vehicle.speed << std::endl;
    int vspeed = (int) vehicle.speed;
    if(Ni.crash_maker) vspeed = 0;

    switch (D) {
        case GoingRight: {
            if (x + vspeed >= msx) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                break;
                //x = -1;
            }
            xy = position_to_cell(msx, msy, x + vspeed, y);
            auto next_cell = ca_matrix.at(xy); //will always be ok.
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            //auto next_vehicle = next_vehicle_iterator != vehicles_map_remote.end()? : &(*next_vehicle_iterator) : nullptr;
            xy = position_to_cell(msx, msy, x+1, y-1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](const auto& v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;

            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case GoingLeft: {

            if (x-vspeed < 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;

                break; x = msx;
            }

            xy = position_to_cell(msx, msy, x-vspeed, y);
            auto next_cell = ca_matrix.at(xy);
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            xy = position_to_cell(msx, msy, x-1, y+1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](auto v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;
            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case GoingUp: {

            if (y-vspeed < 0) { // deletion at boundary, if parallel, ask neighbor
                //vehicles_map_new[y][x] = nullptr;
                break; y = msy;
            }

            xy = position_to_cell(msx, msy, x, y - vspeed);
            auto next_cell = ca_matrix.at(xy);
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            xy = position_to_cell(msx, msy, x-1 == -1 ? msx : x - 1, y - 1 == -1 ? msy : y - 1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](auto v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;
            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case GoingDown: {
            if(y+vspeed >= msy) break;
            xy = position_to_cell(msx, msy, x, y + vspeed);
            auto next_cell = ca_matrix.at(xy);
            const Vehicle* next_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            xy = position_to_cell(msx, msy, x+1 == msx ? 0 : x+1, y+1 == msy ? 0 : y+1);
            auto priority_vehicle = get_or_default(ca_matrix, xy, NoDirection, [](auto v){return v.direction;}) == Rotary ? get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy) : nullptr;
            parallel::in_road_rule(msx, msy, &vehicle, next_cell, next_vehicle, priority_vehicle, ca_matrix, vehicles_map, vehicles_map_remote, vehicles_map_new);

        }
            break;
        case Rotary: { // von neumann neighborhood
            std::pair<int,int> exit, next; std::tie(rid, exit, next) = get_rotary_neighbors(msx, msy, x, y, ca_matrix);
            auto exit_cell  = ca_matrix.at(position_to_cell(msx, msy, exit));
            auto next_cell  = ca_matrix.at(position_to_cell(msx, msy, next));
            xy = position_to_cell(msx,msy, exit);
            const Vehicle* exit_vehicle = get_ptr_vehicle(vehicles_map, vehicles_map_remote, xy);
            if(rid == 3) vehicle.rotary_exit_flag = 1;
            sequential::in_rotary_rule(msx, msy, &vehicle, exit_cell, next_cell, exit_vehicle, vehicles_map_new);
        }
            break;
        case NoDirection:
            std::cerr << "No direction: skipping..." << std::endl;
    }
}
}// end of namespace: parallel


#endif //CA_ROAD_RULES_HPP

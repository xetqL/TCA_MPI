//
// Created by xetql on 20.12.18.
//

#ifndef CA_ROAD_VEHICLE_HPP
#define CA_ROAD_VEHICLE_HPP

#include <tuple>
#include <ostream>
#include "ca_cell.hpp"
#include "communication.hpp"
#define PROBABILITY_EXIT 50
struct Vehicle {
    int gid = 0;
    int lid = 0;
    int rotary_exit_flag;
    std::pair<int, int> position;
    float speed;

    Vehicle(): gid(0), lid(0), rotary_exit_flag(0), position(std::make_pair(0, 0)),  speed(0) {}
    Vehicle(int gid, int lid, int x, int y, float speed) : gid(gid), lid(lid), rotary_exit_flag(rand() % 100 < PROBABILITY_EXIT), position(std::make_pair(x, y)),  speed(speed) {}
    Vehicle(int gid, int lid, std::pair<int, int> pos, float speed) : gid(gid), lid(lid), rotary_exit_flag(rand() % 100 < PROBABILITY_EXIT), position(std::move(pos)),  speed(speed) {}
    Vehicle(int gid, int lid, CA_Cell& c, float speed) : gid(gid), lid(lid), rotary_exit_flag(rand() % 100 < PROBABILITY_EXIT), position(c.position),  speed(speed) {}
    Vehicle(int gid, int lid, const CA_Cell& c, float speed) : gid(gid), lid(lid), rotary_exit_flag(rand() % 100 < PROBABILITY_EXIT), position(c.position),  speed(speed) {}

    friend std::ostream &operator<<(std::ostream &os, const Vehicle &vehicle) {
        int x,y;
        std::tie(x, y) = vehicle.position;
        os << "gid: " << vehicle.gid << " lid: " << vehicle.lid << " rotary_exit_flag: " << vehicle.rotary_exit_flag
           << " position: (" << x << ","<< y <<") speed: " << vehicle.speed;
        return os;
    }

    static CommunicationDatatype register_datatype() {

        MPI_Datatype element_datatype,
                position_datatype,
                oldtype_element[3];
        MPI_Aint offset[3], intex, pos_offset;
        int blockcount_element[3];

        // register particle element type
        MPI_Type_contiguous(2, MPI_INT, &position_datatype); // position
        MPI_Type_commit(&position_datatype);

        blockcount_element[0] = 3; // gid, lid, exit
        blockcount_element[1] = 1; // position <x,y>
        blockcount_element[2] = 1; // speed

        oldtype_element[0] = MPI_INT;
        oldtype_element[1] = position_datatype;
        oldtype_element[2] = MPI_FLOAT;

        MPI_Type_extent(MPI_INT, &intex);
        MPI_Type_extent(position_datatype, &pos_offset);


        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = 3 * intex;
        offset[2] = 3 * intex + pos_offset;

        MPI_Type_struct(3, blockcount_element, offset, oldtype_element, &element_datatype);

        MPI_Type_commit(&element_datatype);

        return CommunicationDatatype(position_datatype, element_datatype);
    }
};


#endif //CA_ROAD_VEHICLE_HPP

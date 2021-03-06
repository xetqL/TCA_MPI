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

    int waiting_time = 0;
    std::pair<int, int> position;
    float speed = 1;
    mutable int rotary_exit_flag;

    Vehicle(): gid(-1), lid(-1), rotary_exit_flag(0), position(std::make_pair(-1, -1)),  speed(1) {}
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

    std::string as_readable_string() {
        return std::to_string(position.first) + ";"+std::to_string(position.second)+";"+std::to_string(waiting_time);
    }

    char as_char() const {
        return '#';
    }

    std::array<double, 2> get_position_as_array() {
        return {(double) position.first, (double) position.second};
    };

    Vehicle drive(int x, int y, float speed) const {
        //usleep((rand() % 500 + 500)); //simulate a more complex Agent's task
        return Vehicle(gid, lid, x, y, speed);
    }

    Vehicle drive(std::pair<int,int> pos, float speed) const {
        int x, y; std::tie(x,y) = pos;
        return this->drive(x,y,speed);
    }

    bool operator==(const Vehicle &rhs) const {
        return position == rhs.position;
    }

    bool operator!=(const Vehicle &rhs) const {
        return !(rhs == *this);
    }

    static CommunicationDatatype register_datatype() {

        MPI_Datatype element_datatype,
                position_datatype,
                oldtype_element[3];
        MPI_Aint offset[3], intex, pos_offset;

        const int number_of_int_elements = 4;
        const int number_of_position_elements = 1;
        const int number_of_speed_elements = 1;

        int blockcount_element[3];

        // register particle element type
        MPI_Type_contiguous(2, MPI_INT, &position_datatype); // position
        MPI_Type_commit(&position_datatype);

        blockcount_element[0] = number_of_int_elements; // gid, lid, exit, waiting_time
        blockcount_element[1] = number_of_position_elements; // position <x,y>
        blockcount_element[2] = number_of_speed_elements; // speed

        oldtype_element[0] = MPI_INT;
        oldtype_element[1] = position_datatype;
        oldtype_element[2] = MPI_FLOAT;

        MPI_Type_extent(MPI_INT, &intex);
        MPI_Type_extent(position_datatype, &pos_offset);


        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = number_of_int_elements * intex;
        offset[2] = number_of_int_elements * intex + number_of_position_elements * pos_offset;

        MPI_Type_struct(3, blockcount_element, offset, oldtype_element, &element_datatype);

        MPI_Type_commit(&element_datatype);

        return CommunicationDatatype(position_datatype, element_datatype);
    }

    bool want_to_exit(int rotary_pos){
        return rotary_exit_flag || rotary_pos == 3;
    }

    void set_exit_flag(int exit) {
        this->rotary_exit_flag = exit;
    }

};


#endif //CA_ROAD_VEHICLE_HPP

//
// Created by xetql on 20.12.18.
//

#ifndef CA_ROAD_CA_CELL_HPP
#define CA_ROAD_CA_CELL_HPP

#include "driving_direction.hpp"

struct CA_Cell {
    DrivingDirection direction;
    std::pair<int,int> position;

    CA_Cell() : direction(DrivingDirection::NoDirection) {}
    CA_Cell(char c) : direction(char_to_driving(c)) {}
    CA_Cell(std::pair<int,int> position, char c) : position(std::move(position)), direction(char_to_driving(c)) {}

};

#endif //CA_ROAD_CA_CELL_HPP

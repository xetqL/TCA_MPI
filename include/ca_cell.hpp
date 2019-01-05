//
// Created by xetql on 20.12.18.
//

#ifndef CA_ROAD_CA_CELL_HPP
#define CA_ROAD_CA_CELL_HPP

#include <ostream>
#include "driving_direction.hpp"

struct CA_Cell {
    DrivingDirection direction = NoDirection;
    std::pair<int,int> position;
    bool source = false;
    mutable int PROBABILITY_GENERATE = rand() % 100;

    CA_Cell() = default;
    CA_Cell(char c) : direction(char_to_driving(c)) {}
    CA_Cell(std::pair<int,int> position, char c) : position(std::move(position)), direction(char_to_driving(c)) {}
    CA_Cell(std::pair<int,int> position, DrivingDirection c) : position(std::move(position)), direction(c) {}

    friend std::ostream &operator<<(std::ostream &os, const CA_Cell &cell) {
        os << "direction: " << driving_to_char(cell.direction) << " position: [" << cell.position.first << ","<<cell.position.second<<"]";
        return os;
    }

    char as_char() const {
        return driving_to_char(direction);
    }

    static CA_Cell& get_fake_cell(){
        static CA_Cell fake_cell = {{-1,-1}, NoDirection};
        return fake_cell;
    }

    bool has_to_generate(int step) const {
        if(step % 10){
            PROBABILITY_GENERATE = rand() % 100;
        }
        return source && rand() % 100 <= PROBABILITY_GENERATE;
    }

};

#endif //CA_ROAD_CA_CELL_HPP

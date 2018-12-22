//
// Created by xetql on 12/22/18.
//

#ifndef CA_ROAD_TCA_IO_HPP
#define CA_ROAD_TCA_IO_HPP

#include <vector>
#include <fstream>
#include <unordered_map>
#include "ca_cell.hpp"
#include "tca_utils.hpp"

std::pair<int, int> read_roadfile(std::string filename, std::vector<std::vector<CA_Cell> > *ca_matrix) {
    std::ifstream f;
    f.open(filename);
    std::string line;
    ca_matrix->clear();
    size_t y = 0;
    while (std::getline(f, line)) {
        size_t x = 0;
        std::vector<char> v(line.begin(), line.end());
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

std::pair<int, int> read_roadfile(std::string filename, std::unordered_map<long long, CA_Cell> *ca_matrix) {
    std::ifstream f;
    f.open(filename);
    std::string line;
    ca_matrix->clear();
    size_t y = 0;
    int msx;
    while (std::getline(f, line)) {
        msx = line.size();
        int x = 0;
        std::vector<char> v(line.begin(), line.end());
        for (char c : v) {
            ca_matrix->operator[](position_to_cell(msx, -1, x, y)) = CA_Cell(std::make_pair(x,y), c);
            x++;
        }
        y++;
    }
    return std::make_pair(msx, y);
}

void print(size_t sx, size_t sy, const std::unordered_map<long long, CA_Cell> &ca_matrix,
           const std::unordered_map<long long, Vehicle> &vehicle_matrix) {
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix.at(position_to_cell(sx, sy, x, y)).direction == NoDirection) std::cout << "*";
            else {
                if (exists(vehicle_matrix, position_to_cell(sx, sy, x, y)))
                    std::cout << "#";
                else std::cout << " ";
            }
        }
        std::cout << std::endl;
    }
}

#endif //CA_ROAD_TCA_IO_HPP

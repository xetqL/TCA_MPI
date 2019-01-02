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
#include "../src/zupply.hpp"

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

void fprint(std::ofstream& out, size_t sx, size_t sy, const std::unordered_map<long long, CA_Cell> &ca_matrix,
           const std::unordered_map<long long, Vehicle> &vehicle_matrix) {
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            auto xy = position_to_cell(sx, sy, x, y);
            if (ca_matrix.at(xy).direction == NoDirection) out << ca_matrix.at(xy).as_char();
            else {
                if (exists(vehicle_matrix, xy))
                    out << vehicle_matrix.at(xy).as_char();
                else out << " ";
            }
        }
        out << std::endl;
    }
}

uint8_t* to_frame(size_t sx, size_t sy, const std::unordered_map<long long, CA_Cell> &ca_matrix,
                          const std::unordered_map<long long, Vehicle> &vehicle_matrix){
    const size_t sz = sy*sx;
    uint8_t *ret = new uint8_t[sz*3];
    for (size_t xy = 0; xy < sz; xy+=3) {
        if (exists(vehicle_matrix, xy)){
            ret[xy]   = 0;
            ret[xy+1] = 0;
            ret[xy+2] = 0;
        }
        else {
            ret[xy]   = 0;
            ret[xy+1] = 0;
            ret[xy+2] = 0;
        }
    }
    /*for (size_t xy = 0; xy < sz; xy++) {
        if (exists(vehicle_matrix, xy))
            ret[xy] = 0;
        else ret[xy] = 0;
    }*/
    return ret;
}

zz::Image zzframe(std::string fname, size_t sx, size_t sy, const std::unordered_map<long long, CA_Cell> &ca_matrix,
                  const std::unordered_map<long long, Vehicle> &vehicle_matrix){
    const size_t sz = sy*sx;
    zz::Image ret(sy, sx, 3);
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            auto xy = position_to_cell(sx, sy, x, y);
            if (exists(vehicle_matrix, xy)) {
                ret(y, x, 0) = 0;
                ret(y, x, 1) = 0;
                ret(y, x, 2) = 0;
            } else {
                ret(y, x, 0) = 255;
                ret(y, x, 1) = 255;
                ret(y, x, 2) = 255;
            }
        }
    }
    /*for (size_t xy = 0; xy < sz; xy++) {
        if (exists(vehicle_matrix, xy))
            ret[xy] = 0;
        else ret[xy] = 0;
    }*/
    return ret;
}


#endif //CA_ROAD_TCA_IO_HPP

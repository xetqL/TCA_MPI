#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <zconf.h>
#include <memory>

#define SPATIAL_DISCRETISATION_X = 7.5 //the average length a conventional vehicle occupies in a closely packed jam (and as such, its width is neglected),
#define TEMPORAL_DISCRETISATION = 1.0 //typical driver’s reaction time


enum CA_State { Empty, Car, Wall }; // with int: 0, 1, 2
enum DrivingDirection { GoingLeft, GoingRight, GoingUp, GoingDown, Rotary, NoDirection }; // with int: 0, 1, 2, 3, 4

using namespace std;

int char_to_state(char c) {
    switch (c){
        case 'W': return 2;
        case 'C': return 1;
        default : return 0;
    }
}
DrivingDirection char_to_driving(char c){
    switch (c) {
        case 'L': return GoingLeft;
        case 'R': return GoingRight;
        case 'U': return GoingUp;
        case 'D': return GoingDown;
        default : return NoDirection;
    }
}
struct DrivingElement {
    int rotary_exit_flag;
    DrivingElement(int flag) : rotary_exit_flag(flag) {

    }
};
struct CA_Cell {
    std::tuple<int, int> position;
    //std::unique_ptr<DrivingElement> element;
    int state;
    DrivingDirection direction;
    CA_Cell(): position(std::make_tuple(0, 0)), state(2), direction(DrivingDirection::NoDirection) {  }
    CA_Cell(int x, int y, char c) : position(std::make_tuple(x, y)), state(char_to_state(c)), direction(char_to_driving(c))
    {}
};

double speed(int delta_x_meters, int delta_y_seconds){
    return (delta_x_meters / 1000.0) / (delta_y_seconds / 3600.0);
}

void init_empty_ca_matrix(vector<vector<CA_Cell> > *ca_matrix) {
    size_t sx = ca_matrix->size(), sy = ca_matrix->operator[](0).size();
    for(size_t x = 0; x < sx; x++) {
        for(size_t y = 0; y < sy; y++) {
            ca_matrix->operator[](x)[y].position = std::make_tuple(x, y);
        }
    }
}

template<class A>
void update(size_t sx, size_t sy, A &rule, vector<vector<CA_Cell> > *ca_matrix){
    vector<vector<CA_Cell> > cp_ca_matrix = *ca_matrix; //deep copy
    for(size_t y = 0; y < sy; ++y) {
        for(size_t x = 0; x < sx; ++x) {
            ca_matrix->operator[](y)[x].state = rule.apply(x, y, cp_ca_matrix);
        }
    }
}

class Rule184 {
    int msx, msy;
public:
    Rule184(int sx, int sy) : msx(sx), msy(sy){}

    int apply(int x, int y, const vector<vector<CA_Cell> > &ca_matrix) {
        //int x, y;
        //std::tie(x, y) = position;
        CA_Cell Ni = ca_matrix.at(y).at(x);
        if(Ni.state == 2) return 2; // don't treat them
        auto D = Ni.direction;

        /** Neighborhood **/
        int prevNi_state, nextNi_state;
        switch (D) {
            case GoingRight:
                if (x+1 < msx){  nextNi_state = ca_matrix.at(y).at(x+1).state; } else { nextNi_state = 0; }
                if (x-1 >= 0) {  prevNi_state = ca_matrix.at(y).at(x-1).state; } else { prevNi_state = rand() % 2; }
                break;
            case GoingLeft:
                if (x+1 < msx){ prevNi_state = ca_matrix.at(y).at(x+1).state; } else { prevNi_state = 1; }
                if (x-1 >= 0) { nextNi_state = ca_matrix.at(y).at(x-1).state; } else { nextNi_state = 0; }
                break;
            case GoingUp:
                if (y-1 >= 0){ prevNi_state = ca_matrix.at(y-1).at(x).state; }  else { prevNi_state = 1; }
                if (y+1 < msy){ nextNi_state = ca_matrix.at(y+1).at(x).state; } else { nextNi_state = 0; }
                break;
            case GoingDown:
                if (y-1 >= 0) { nextNi_state = ca_matrix.at(y-1).at(x).state; } else { nextNi_state = 0; }
                if (y+1 < msy){ prevNi_state = ca_matrix.at(y+1).at(x).state; } else { prevNi_state = 1; }
                break;
            /*case Rotary: // von neumann neighborhood
                auto north = ca_matrix.at(y-1).at(x); auto south = ca_matrix.at(y+1).at(x);
                auto east  = ca_matrix.at(y).at(x+1); auto west  = ca_matrix.at(y-1).at(x-1);
                break;*/
            default:
                throw std::runtime_error("unknown direction.");
        }
        return prevNi_state * (1 - Ni.state) + Ni.state * nextNi_state;
    }
};

std::tuple<int,int> read_roadfile(std::string filename, vector<vector<CA_Cell> > *ca_matrix){
    std::ifstream f;
    f.open(filename);
    std::string line;
    ca_matrix->clear();
    size_t y = 0;
    while ( std::getline(f, line) ) {
        size_t x = 0;
        vector<char> v(line.begin(), line.end());
        std::vector<CA_Cell> road;
        for(char c : v) {
            road.emplace_back(x, y, c);
            x++;
        }
        ca_matrix->push_back(road);
        y++;
    }
    return std::make_tuple(ca_matrix->operator[](0).size(), y);
}

void randomize_cars_position(size_t sx, size_t sy, vector<vector<CA_Cell> > *ca_matrix){
    for(size_t x = 0; x < sx; x++) {
        for (size_t y = 0; y < sy; y++) {
            if(ca_matrix->operator[](y)[x].state == 0)
                if(rand() % 100 > 55)
                    ca_matrix->operator[](y)[x].state = 1;
                    //ca_matrix->operator[](y)[x].element = std::make_unique<DrivingElement>(1);
        }
    }
}

void print_CA(const vector<vector<CA_Cell> >& ca_matrix) {
    for(auto l : ca_matrix){
        for(auto el : l){
            switch (el.state){
                case 0:
                    switch(el.direction) {
                        case GoingRight: std::cout << "R"; break;
                        case GoingDown: std::cout  << "D"; break;
                        case GoingLeft: std::cout  << "L"; break;
                        case GoingUp: std::cout  << "U"; break;
                    }
                    break;
                case 1: std::cout << 'c'; break;
                case 2: std::cout << 'W'; break;
            }
        }
        std::cout << std::endl;
    }
}

int main(int argc, char** argv) {
    int SIZE_X,SIZE_Y;
    vector<vector<CA_Cell> >  ca_matrix;

    std::tie(SIZE_X, SIZE_Y) = read_roadfile("roads.txt", &ca_matrix);
    randomize_cars_position(SIZE_X, SIZE_Y, &ca_matrix);
    Rule184 r(SIZE_X, SIZE_Y);

    int step = 0; print_CA(ca_matrix);
    while(step < 100) {
        system("clear");
        update(SIZE_X, SIZE_Y, r, &ca_matrix);
        print_CA(ca_matrix);
        sleep(2); step++;
    }

    //print_CA(ca_matrix);
    return 0;
}
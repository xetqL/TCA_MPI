//
// Created by xetql on 12/21/18.
//

#ifndef CA_ROAD_TCA_UTILS_HPP
#define CA_ROAD_TCA_UTILS_HPP

#include <map>
#include <unordered_map>
#include <algorithm>

template<class key, class val, class retval, class RetFunc>
inline retval get_or_default(const std::unordered_map<key, val>& map, const key& k1, const retval& def,  RetFunc f){
    if(map.find(k1) == map.cend())
        return def;
    return f(map.at(k1));
}

template<class key, class val, class retval>
inline retval get_or_default(const std::unordered_map<key, val>& map, const key& k1, const retval& def) {
    if(map.find(k1) == map.cend())
        return def;
    return map.at(k1);
}

template<class key, class val, class retval>
inline retval get_or_default(const std::unordered_map<key, val> map, const key k1, const retval def) {
    if(map.find(k1) == map.cend())
        return def;
    return map.at(k1);
}

template<class Map, class key>
inline bool exists(const Map &map, const key k1){
    return map.find(k1) != map.cend();
}

inline double distance2(const Vehicle& v1, const Vehicle& v2){
    const double dx = v1.position.first-v2.position.first;
    const double dy = v1.position.second-v2.position.second;
    return dx*dx+dy*dy;
}

inline const Vehicle* get_ptr_vehicle(const std::unordered_map<long long, Vehicle>& vehicles_map,
                         const std::unordered_map<long long, Vehicle>& vehicles_map_remote,
                         long long xy) {

    if (exists(vehicles_map, xy)){
        return &vehicles_map.at(xy);
    } else if(exists(vehicles_map_remote, xy)) {
        return &vehicles_map_remote.at(xy);
    } else {
        return nullptr;
    }
}

template<class Map, class key>
inline int exists(const Map &map1, const Map &map2, const key k1){
    if(map1.find(k1) != map1.cend())
        return 1;
    if(map2.find(k1) != map2.cend())
        return 2;

    return 0;
}

inline long long position_to_cell(int msx, int msy, const std::pair<int, int> & position) {
    return position.first + msx * position.second;
}

inline long long position_to_cell(int msx, int msy, const int x, const int y) {
    return x + msx * y;
}

inline std::pair<int, int> cell_to_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msx);
}


void randomize_cars_position(size_t sx, size_t sy, const std::vector<std::vector<CA_Cell> > &ca_matrix, std::vector<Vehicle> *vehicles,
                             std::vector<std::vector<Vehicle *>> *vehicle_matrix) {
    int gid = 0;
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix[y][x].direction != NoDirection) {
                if (rand() % 100 > 66) {
                    vehicles->emplace_back(gid, gid, x, y, 1);
                    vehicle_matrix->operator[](y)[x] = &(*(vehicles->end() - 1));
                    gid++;
                }
            }
        }
    }
}

void randomize_cars_position(size_t sx, size_t sy, const std::unordered_map<long long, CA_Cell> &ca_matrix,
                             std::unordered_map<long long, Vehicle>& vehicle_matrix) {

    int gid = 0;
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix.at(position_to_cell(sx, sy, x, y)).direction != NoDirection) {
                if (rand() % 100 > 50) {
                    vehicle_matrix[position_to_cell(sx, sy, x, y)] = Vehicle(gid, gid, x, y, 1);
                    gid++;
                }
            }
        }
    }
}

void randomize_cars_position(int nb_car, size_t sx, size_t sy, const std::vector<std::vector<CA_Cell> > &ca_matrix, std::vector<Vehicle> *vehicles,
                             std::vector<std::vector<Vehicle *>> *vehicle_matrix) {
    int gid = 0;
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix[y][x].direction != NoDirection) {
                if (rand() % 100 > 10 && gid < nb_car) {
                    vehicles->emplace_back(gid, gid, x, y, 1);
                    vehicle_matrix->operator[](y)[x] = &(*(vehicles->end() - 1));
                    gid++;
                }
            }
        }
    }
}

void print(size_t sx, size_t sy, const std::vector<std::vector<CA_Cell> > &ca_matrix,
           const std::vector<std::vector<Vehicle *>> &vehicle_matrix) {
    for (size_t y = 0; y < sy; y++) {
        for (size_t x = 0; x < sx; x++) {
            if (ca_matrix.at(y).at(x).direction == NoDirection) std::cout << "*";
            else {
                if (vehicle_matrix.at(y).at(x) == nullptr)
                    std::cout << " ";
                else std::cout << "#";
            }
        }
        std::cout << std::endl;
    }
}

int count_car(int msx, int msy, std::vector<std::vector<Vehicle *> >& vehicle_matrix){
    int c = 0;
    for(int i = 0; i < msy; i++){
        for(int j = 0; j < msx; j++){
            if(vehicle_matrix[i][j] != nullptr) c++;
        }
    }
    return c;
}

int count_car(int msx, int msy, const std::unordered_map<long long, Vehicle> &vehicle_matrix){
    return vehicle_matrix.size();
}

int char_to_state(char c) {
    switch (c) {
        case 'W':
            return 2;
        case 'C':
            return 1;
        default :
            return 0;
    }
}

std::vector<Vehicle> to_vec(std::unordered_map<long long, Vehicle> &vehicle_matrix){
    std::vector<Vehicle> vec;
    vec.reserve(vehicle_matrix.size());
    size_t i = 0;
    for(auto& v : vehicle_matrix){
        vec.push_back(std::move(v.second));
        i++;
    }
    return vec;
}

std::unordered_map<long long, Vehicle> to_map(int msx, int msy, std::vector<Vehicle>& vec) {
    std::unordered_map<long long, Vehicle> vehicle_matrix;
    vehicle_matrix.reserve(vec.size());
    size_t vsize = vec.size();
    long long xy;
    for(size_t i = 0; i < vsize; i++){
        xy = position_to_cell(msx, msy, vec[i].position);
        vehicle_matrix[xy] = std::move(vec[i]);
    }
    return vehicle_matrix;
}

std::tuple<std::unordered_map<long long, CA_Cell>, std::vector<int>, std::vector<int>>  generate_random_manhattan (const size_t SIZE_X, const size_t SIZE_Y) {
    const long long MANHATTAN_LENGTH = SIZE_X * SIZE_Y;
    std::unordered_map<long long, CA_Cell> manhattan;
    manhattan.reserve(MANHATTAN_LENGTH);

    std::vector<int> roads_position_x;
    for(int last = (rand() % 7) + 1; (last+2) < SIZE_X; last += (rand() % 7) + 3){
        roads_position_x.push_back(last);
    }

    std::vector<int> roads_position_y;
    for(int last = (rand() % 7) + 1; (last + 2) < SIZE_Y; last += (rand() % 7) + 3) {
        roads_position_y.push_back(last);
    }

    for(long long i = 0; i < MANHATTAN_LENGTH; ++i) {
        int x, y; std::tie(x,y) = cell_to_position(SIZE_X, SIZE_Y, i);
        manhattan[i] = { std::make_pair(x, y), NoDirection };
    }

    for(const int road_x : roads_position_x){
        for(int road_y = 0; road_y < SIZE_Y; road_y++) {
            long long cell_idx = position_to_cell(SIZE_X, SIZE_Y, std::make_pair(road_x, road_y));
            manhattan[cell_idx].direction = GoingDown;
            manhattan[cell_idx+1].direction = GoingUp;
        }
    }

    for(const int road_y : roads_position_y){
        for(int road_x = 0; road_x < SIZE_X; road_x++) {
            long long cell_idx = position_to_cell(SIZE_X, SIZE_Y, std::make_pair(road_x, road_y));
            if(manhattan[cell_idx].direction != NoDirection)
                manhattan[cell_idx].direction = Rotary;
            else
                manhattan[cell_idx].direction = GoingLeft;

            if(manhattan[cell_idx+SIZE_X].direction != NoDirection)
                manhattan[cell_idx+SIZE_X].direction = Rotary;
            else
                manhattan[cell_idx+SIZE_X].direction = GoingRight;
        }
    }
    return std::make_tuple(manhattan, roads_position_x, roads_position_y);
}

void create_random_left_sources(int n, const size_t SIZE_X, const size_t SIZE_Y, std::vector<int> left, std::unordered_map<long long, CA_Cell> *ca_matrix) {
    const int X = 0;
    std::random_shuffle(left.begin(), left.end());
    n = n > left.size() ? left.size() : n;
    for(; n > 0; n--) {
        const int Y = left[n-1]+1;
        const long long xy = position_to_cell(SIZE_X, SIZE_Y, X, Y);
        ca_matrix->at(xy).source = true;
    }
}

std::tuple<int, int, int, int> get_geometry_from_vehicles(
        int rank, Zoltan_Struct* zz, const std::unordered_map<long long, Vehicle>& my_vehicles, const size_t SIZE_X, const size_t SIZE_Y) {
    int x, y, minx= SIZE_X, miny = SIZE_Y, maxx=-1,maxy=-1;
    //create boundaries from vehicles
    for(const auto& v : my_vehicles) {
        std::tie(x, y) = cell_to_position(SIZE_X, SIZE_Y, v.first);
        minx = std::min(x, minx); miny = std::min(y, miny);
        maxx = std::max(x, maxx); maxy = std::max(y, maxy);
    }
    //then grow from those boundaries to the edges of my own area
    //using linear search (but maybe, binary search could provide better perf.) TODO: test with binary search
    int PE, part;
    std::array<double, 2> pos_in_double;
    for(int i = minx; i > 0; --i) {
        pos_in_double = {(double) i, (double) miny};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            minx = i;
    }

    for(int i = miny; i > 0; --i) {
        pos_in_double = {(double) minx, (double) i};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            minx = i;
    }

    for(int i = maxx; i < SIZE_X; ++i) {
        pos_in_double = {(double) i, (double) maxy};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            maxx = i;
    }

    for(int i = maxy; i < SIZE_Y; ++i) {
        pos_in_double = {(double) maxx, (double) i};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            maxx = i;
    }

    return std::make_tuple(minx, maxx, miny, maxy);
}
std::tuple<int, int, int, int> get_geometry_from_vehicles(
                    int rank, Zoltan_Struct* zz,
                    const std::vector<Vehicle>& my_vehicles,
                    const size_t SIZE_X, const size_t SIZE_Y) {
    int x, y, minx= SIZE_X, miny = SIZE_Y, maxx=-1,maxy=-1;
    //create boundaries from vehicles
    for(const auto& v : my_vehicles) {
        std::tie(x, y) = v.position;
        minx = std::min(x, minx); miny = std::min(y, miny);
        maxx = std::max(x, maxx); maxy = std::max(y, maxy);
    }
    //then grow from those boundaries to the edges of my own area
    //using linear search (but maybe, binary search could provide better perf.) TODO: test with binary search
    int PE, part;
    std::array<double, 2> pos_in_double;
    for(int i = minx; i > 0; --i) {
        pos_in_double = {(double) i, (double) miny};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            minx = i;
    }

    for(int i = miny; i > 0; --i) {
        pos_in_double = {(double) minx, (double) i};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            minx = i;
    }

    for(int i = maxx; i < SIZE_X; ++i) {
        pos_in_double = {(double) i, (double) maxy};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            maxx = i;
    }

    for(int i = maxy; i < SIZE_Y; ++i) {
        pos_in_double = {(double) maxx, (double) i};
        Zoltan_LB_Point_PP_Assign(zz, &pos_in_double.front(), &PE, &part);
        if(PE == rank) // my point
            maxx = i;
    }

    return std::make_tuple(minx, maxx, miny, maxy);
}

inline bool inside_the_borders(const Vehicle& vehicle, const std::tuple<int, int, int, int>& my_area, const int border_size=1) {
    return (vehicle.position.first  <= std::get<0>(my_area)+border_size)|| (std::get<1>(my_area)-border_size <= vehicle.position.first) ||
           (vehicle.position.second <= std::get<2>(my_area)+border_size)|| (std::get<3>(my_area)-border_size <= vehicle.position.second);
}


void set_crash_maker_on_out_roads(bool crash, const size_t SIZE_X, const size_t SIZE_Y, const std::tuple<int, int, int, int>& geom,
        std::unordered_map<long long, CA_Cell>* _ca_matrix){
    std::unordered_map<long long, CA_Cell>& ca_matrix = *(_ca_matrix);
    int minx,maxx, miny,maxy; std::tie(minx,maxx, miny,maxy) = geom;
    int X, Y;

    Y=miny+1;
    for(X = minx; X < maxx; ++X) {
        const long long xy = position_to_cell(SIZE_X, SIZE_Y, X, Y);
        auto& cell = ca_matrix.at(xy);
        if(cell.direction == GoingUp)
            cell.crash_maker = crash;
        else if (cell.direction == Rotary) {
            int tmp_Y = Y;
            CA_Cell *tmp_cell;
            do {
                tmp_Y++;
                const long long xy = position_to_cell(SIZE_X, SIZE_Y, X, tmp_Y);
                tmp_cell = &ca_matrix.at(xy);

            } while(tmp_cell->direction == Rotary);
            if(tmp_cell->direction == GoingUp) tmp_cell->crash_maker = crash;
        }
    }
    Y=maxy-1;
    for(X = minx; X < maxx; ++X) {
        const long long xy = position_to_cell(SIZE_X, SIZE_Y, X, Y);
        auto& cell = ca_matrix.at(xy);
        if(cell.direction == GoingDown)
            cell.crash_maker = crash;
        else if (cell.direction == Rotary) {
            int tmp_Y = Y;
            CA_Cell *tmp_cell;
            do {
                tmp_Y--;

                const long long xy = position_to_cell(SIZE_X, SIZE_Y, X, tmp_Y);
                tmp_cell = &ca_matrix.at(xy);
            } while(tmp_cell->direction == Rotary);
            if(tmp_cell->direction == GoingDown)  tmp_cell->crash_maker = crash;
        }
    }

    X=minx+1;
    for(Y = miny; Y < maxy; ++Y) {
        const long long xy = position_to_cell(SIZE_X, SIZE_Y, X, Y);
        auto& cell = ca_matrix.at(xy);
        if(cell.direction == GoingLeft)
            cell.crash_maker = crash;
        else if (cell.direction == Rotary) {
            int tmp_X = X;
            CA_Cell *tmp_cell;
            do {
                tmp_X++;
                const long long xy = position_to_cell(SIZE_X, SIZE_Y, tmp_X, Y);
                tmp_cell = &ca_matrix.at(xy);
            } while(tmp_cell->direction == Rotary);

            if(tmp_cell->direction == GoingLeft) tmp_cell->crash_maker = crash;
        }
    }

    X=maxx-1;
    for(Y = miny; Y < maxy; ++Y) {
        const long long xy = position_to_cell(SIZE_X, SIZE_Y, X, Y);
        auto& cell = ca_matrix.at(xy);
        if(cell.direction == GoingRight)
            cell.crash_maker = crash;
        else if (cell.direction == Rotary) {
            int tmp_X = X;
            CA_Cell *tmp_cell;
            do {
                tmp_X--;
                const long long xy = position_to_cell(SIZE_X, SIZE_Y, tmp_X, Y);
                tmp_cell = &ca_matrix.at(xy);
            } while(tmp_cell->direction == Rotary);
            if(tmp_cell->direction == GoingRight) tmp_cell->crash_maker = crash;
        }
    }
}
#endif //CA_ROAD_TCA_UTILS_HPP

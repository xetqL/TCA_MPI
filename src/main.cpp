#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <zconf.h>
#include <memory>
#include <map>
#include <cassert>
#include <algorithm>

#include "../include/rules.hpp"
#include "../include/zoltan_fn.hpp"
#include "../include/unloading_model.hpp"
#include "../include/tca_io.hpp"
#include "../include/gif.h"

#include "zupply.hpp"
#include "../include/window.hpp"
#include "../include/metric.hpp"

#define SPATIAL_DISCRETISATION_X = 7.5 //the average length a conventional vehicle occupies in a closely packed jam (and as such, its width is neglected),
#define TEMPORAL_DISCRETISATION = 1.0 //typical driver’s reaction time

#define START_TIMING(v) \
auto v = MPI_Wtime()

#define RESTART_TIMING(v) \
v = MPI_Wtime() - v

#define STOP_TIMING(v) \
v = MPI_Wtime() - v

#define PAR_START_TIMING(v, comm) \
MPI_Barrier(comm); \
auto v = MPI_Wtime()

#define PAR_RESTART_TIMING(v, comm) \
std::cout <<"Intermediate result " << v << std::endl;\
MPI_Barrier(comm); \
v = MPI_Wtime() - v

#define PAR_STOP_TIMING(v, comm) \
MPI_Barrier(comm); \
v = MPI_Wtime() - v

using namespace std;
using namespace tca;


#define UNLOADING_MODEL 1 


unordered_map<long long, Vehicle> update(const int msx, const int msy,
                                         const unordered_map<long long,  CA_Cell> &ca_matrix,
                                         const unordered_map<long long,  Vehicle> &vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new;
    for (const auto &v : vehicles_map)
        sequential::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, vehicles_map_new);

    return vehicles_map_new;
}

unordered_map<long long, Vehicle> parallel_update_new_model(const int msx, const int msy,
                                                  const unordered_map<long long,  CA_Cell> &ca_matrix,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map_top,
                                                  const unordered_map<long long,  Vehicle> &remote_vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new, remote_vehicles_map_merged;
    remote_vehicles_map_merged.insert(remote_vehicles_map.begin(), remote_vehicles_map.end());
    remote_vehicles_map_merged.insert(vehicles_map_top.begin(),    vehicles_map_top.end());
    for (const auto &v : vehicles_map)
        parallel::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, remote_vehicles_map_merged, &vehicles_map_new);

    return vehicles_map_new;
}

unordered_map<long long, Vehicle> parallel_update(const int msx, const int msy,
                                                  const unordered_map<long long,  CA_Cell> &ca_matrix,
                                                  const unordered_map<long long,  Vehicle> &vehicles_map,
                                                  const unordered_map<long long,  Vehicle> &remote_vehicles_map) {
    unordered_map<long long, Vehicle> vehicles_map_new;
    for (const auto &v : vehicles_map)
        parallel::apply_rule184(msx, msy, ca_matrix, v.second, vehicles_map, remote_vehicles_map, &vehicles_map_new);

    return vehicles_map_new;
}
std::vector<const CA_Cell*> get_my_cells(const int msx, const int msy, const int my_rank, Zoltan_Struct* zz, const unordered_map<long long,  CA_Cell> &ca_matrix) {
    int _x,_y, PE;
    double coords[2];
    std::vector<const CA_Cell*> res;
    for(auto& c : ca_matrix){
        std::tie(_x,_y) = cell_to_position(msx, msy, c.first);
        coords[0] = _x; coords[1] = _y;
        Zoltan_LB_Point_Assign(zz, coords, &PE);
        if(PE==my_rank) {
            res.push_back(&(c.second));
        }
    }
    return res;
}
template<class Predicate>
std::vector<const CA_Cell*> get_my_cells(const int msx, const int msy, const int my_rank, Zoltan_Struct* zz, const unordered_map<long long,  CA_Cell> &ca_matrix, Predicate pred) {
    int _x,_y, PE;
    double coords[2];
    std::vector<const CA_Cell*> res;
    for(auto& c : ca_matrix){
        if(pred(c.second)){
            std::tie(_x,_y) = cell_to_position(msx, msy, c.first);
            coords[0] = _x; coords[1] = _y;
            Zoltan_LB_Point_Assign(zz, coords, &PE);
            if(PE==my_rank) {
                res.push_back(&(c.second));
            }
        }
    }
    return res;
}

template<typename Iter, typename RandomGenerator>
        Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}
template<typename T>
T select_randomly(T start, T end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

void generate_vehicles(const int step, const int my_rank,
                       const int msx, const int msy, Zoltan_Struct* zz,
                       const unordered_map<long long,  CA_Cell> &ca_matrix,
                             unordered_map<long long,  Vehicle> *_vehicles_map, int random_gen_cnt = 5) {
    int wsize; MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    auto& vehicles_map = *(_vehicles_map);

    auto my_cells = get_my_cells(msx, msy, my_rank, zz, ca_matrix, [&vehicles_map, &msx, &msy](auto cell){
        auto pos = position_to_cell(msx,msy, cell.position);
        return !exists(vehicles_map, pos) && cell.direction != Rotary;
    });

    std::vector<std::pair<long long, const CA_Cell*>> cells;

    for(long long i = 0; i < msx; i+=1) {
        if(ca_matrix.at(i + 0*msx).source) cells.emplace_back(i + 0*msx, &ca_matrix.at(i + 0*msx));
        if(ca_matrix.at(i + (msy-1)*msx).source) cells.emplace_back(i + (msy-1)*msx, &ca_matrix.at(i + (msy-1)*msx));
    }
    for(long long i = 0; i < msy; i+=1) {
        if (ca_matrix.at(i * msx).source) cells.emplace_back(i * msx, &ca_matrix.at(i * msx));
        if (ca_matrix.at((msx - 1) + i * msx).source) cells.emplace_back((msx - 1) + i * msx, &ca_matrix.at((msx - 1) + i * msx));
    }
    //std::copy_if(ca_matrix.cbegin(), ca_matrix.cend(), std::back_inserter(cells), [](auto cell){return cell.second.source;});
    //add new vehicles
    std::array<double, 2> pos = {0, 0};
    int PE, id = 0;
    for (auto &cell : cells) {
        std::tie(pos[0], pos[1]) = cell_to_position(msx, msy, cell.first);
        Zoltan_LB_Point_Assign(zz, &pos.front(), &PE);
        if(my_rank == PE && !exists(vehicles_map, cell.first) && cell.second->has_to_generate(step)) {
            vehicles_map[cell.first] = Vehicle(msx * msy + step * (cells.size()*wsize)+id+my_rank*cells.size(), vehicles_map.size() + 1, (int) pos[0], (int) pos[1], 1);
        }
    }

    int nb_rand_cells = random_gen_cnt > (int) my_cells.size() ? (int) my_cells.size() : random_gen_cnt;
    for(int i = 0; i < nb_rand_cells; ++i) {
        auto& cell = *select_randomly(my_cells.begin(), my_cells.end());
        auto cidx  = position_to_cell(msx,msy,cell->position);
        if(cell->direction != Rotary &&  cell->direction != NoDirection && !exists(vehicles_map, cidx) )
            vehicles_map[cidx] = Vehicle(msx*msy+step*(cells.size()*wsize)+id+my_rank*cells.size()+100, vehicles_map.size() + 1, (int) cell->position.first, cell->position.second, 1);
    }

    for(int i = nb_rand_cells; i < 0 && !vehicles_map.empty(); ++i) {
        auto cell = select_randomly(vehicles_map.begin(), vehicles_map.end());
        vehicles_map.erase(cell);
    }
}

int main(int argc, char **argv) {
    int wsize, rank;
    srand(1);
    MPI_Init(&argc, &argv);
    MPI_Comm bottom = MPI_COMM_WORLD;
    MPI_Comm_size(bottom, &wsize);
    MPI_Comm_rank(bottom, &rank);
    int SIZE_X, SIZE_Y;
    std::ofstream out;

    unordered_map<long long, CA_Cell> ca_matrix;
    const std::string prefix_fname(argv[1]);
    const int MAX_STEP = std::atoi(argv[2]);
    std::vector<int> road_y_pos;

    SIZE_X = std::atoi(argv[3]);
    SIZE_Y = std::atoi(argv[4]);
    const int MAX_VC = SIZE_X*SIZE_Y;
    std::tie(ca_matrix, std::ignore, road_y_pos) = generate_random_manhattan(SIZE_X, SIZE_Y);

    if(!rank) std::cout << "End of map generation\nStarting computations..." << std::endl;

    unordered_map<long long, Vehicle> vehicle_matrix;

    auto datatype = Vehicle::register_datatype();
    if(!rank) {
        randomize_cars_position(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
    }
    auto vehicles = to_vec(vehicle_matrix);
    std::cout << vehicles.size() << std::endl;

    if(!rank){
        auto img = zzframe( SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix);
        img.save("out0.jpg");
    }
    auto zz = zoltan_create_wrapper(ENABLE_AUTOMATIC_MIGRATION, MPI_COMM_WORLD);

    zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);

    auto geom = get_geometry_from_vehicles(rank, zz, vehicles, SIZE_X, SIZE_Y); //wont work with HSFC
    create_random_sources(rank, -1, SIZE_X, SIZE_Y, geom, &ca_matrix); //wont work with HSFC

    std::vector<int> crash(0);
    if(!rank) {
        crash = set_crash_maker_on_out_roads(true, SIZE_X, SIZE_Y, geom, &ca_matrix);
    }

    //sharing all crash positions alltoall
    update_crash_positions(true, SIZE_X, SIZE_Y, crash, &ca_matrix, bottom);

    std::vector<Vehicle> top_vehicles;

#if UNLOADING_MODEL > 0
    if(!rank) std::cout << "UNLOADING MODEL APPLIED" << std::endl;
    double dydx = 0;
    auto model_state = load_balancing::esoteric::init_unloading_model(0, dydx, rank, vehicles, bottom);
    if(model_state.state != load_balancing::esoteric::MODEL_STATE::init) {
        std::cout << "no esoteric call" << std::endl;
    } else std::cout << rank << " Next call in " << model_state.sigma << "iterations "<<std::endl;
    auto zztop = load_balancing::esoteric::start_unloading_model(&vehicles , &top_vehicles, model_state, datatype.elements_datatype, bottom);
    //load_balancing::esoteric::UnloadingModelLocalState model_state;
#endif
    int step = 0;

    SlidingWindow<int> window(5);
    constexpr std::array<int, 5> it = {1,2,3,4,5};
    int rcnt_array[2] = {6, -2}, rcnt = 6, rcnti = -1;
    bool model_stopped = true;
    while (step < MAX_STEP) {
        if(step % 100 == 0) rcnti = (rcnti+1) % 2;
        rcnt = rcnt_array[rcnti];
        window.add(vehicles.size());
        auto dydx = get_slope<float>(it, window.data_container);
        //std::cout << rank << " : " << dydx <<  std::endl;
        if(!rank) {
            //std::for_each(window.data_container.cbegin(), window.data_container.cend(), [](auto v){ std::cout << v << " "; });
            //std::cout << std::endl;
        }
#if UNLOADING_MODEL > 0
        //std::cout << model_stopped << std::endl;
        if(model_stopped) {
            model_state = load_balancing::esoteric::init_unloading_model(step, dydx, rank, vehicles, bottom);
            if(model_state.state == load_balancing::esoteric::init) {
                zztop = load_balancing::esoteric::start_unloading_model(&vehicles , &top_vehicles, model_state, datatype.elements_datatype, bottom);
                if(!rank)
                    std::cout << "Restarting model\n-> Next call in " << model_state.sigma << "iterations "<<std::endl;
                model_stopped = false;
            } else {
                std::cout << "no esoteric call" << std::endl;
                std::cout << vehicles.size() << std::endl;
            }
        }
#endif
        MPI_Barrier(bottom);
        PAR_START_TIMING(step_time, bottom);
        int recv, sent;
        /*************************************Start parallel exchange********************************************/
        PAR_START_TIMING(comm_time, bottom);
#if UNLOADING_MODEL > 0
        std::vector<Vehicle> remote_data;

        remote_data = load_balancing::esoteric::exchange(zz, zztop, &vehicles, &top_vehicles, &recv, &sent, geom, model_state.increasing_cpus, datatype.elements_datatype, bottom, 1.0);
#else
        auto remote_data =  zoltan_exchange_data(zz, &vehicles, &recv, &sent, datatype.elements_datatype, bottom,  1.2);
#endif
        PAR_STOP_TIMING(comm_time, bottom);
        // Stop parallel exchange

        /*************************************Start parallel computation*****************************************/
        PAR_START_TIMING(computation_time, bottom);
        auto vehicle_matrix_remote = to_map(SIZE_X, SIZE_Y, remote_data);
        auto vehicle_matrix_bottom = to_map(SIZE_X, SIZE_Y, vehicles);
        auto vehicle_matrix_top    = to_map(SIZE_X, SIZE_Y, top_vehicles);

        auto updated_vehicle_matrix1 = parallel_update_new_model(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_bottom, vehicle_matrix_top, vehicle_matrix_remote);
        auto updated_vehicle_matrix2 = parallel_update_new_model(SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_top, vehicle_matrix_bottom, vehicle_matrix_remote);

        if(!rank) generate_vehicles(step, rank, SIZE_X, SIZE_Y, zz, ca_matrix, &updated_vehicle_matrix1, rcnt);

        vehicles     = to_vec(updated_vehicle_matrix1);
        top_vehicles = to_vec(updated_vehicle_matrix2);

        PAR_STOP_TIMING(computation_time, bottom);

        // Stop parallel computation
        /********************************Start load balancing and migration**************************************/
        PAR_START_TIMING(migrate_time, bottom);
#if UNLOADING_MODEL > 0

        load_balancing::esoteric::migrate(zz, zztop, &vehicles, &top_vehicles, model_state.increasing_cpus, datatype.elements_datatype, bottom);

#else
        zoltan_migrate_particles(zz, &vehicles, datatype.elements_datatype, bottom);
        //zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);
#endif
        PAR_STOP_TIMING(migrate_time, bottom);
        // Stop load balancing and migration
#if UNLOADING_MODEL > 0
        if(!model_stopped) {
            model_stopped = load_balancing::esoteric::stop_unloading_model(step, zz, zztop, &vehicles, &top_vehicles, model_state, bottom);
            if(model_stopped) {
                zztop = nullptr;
                zoltan_load_balance(&vehicles, zz, ENABLE_AUTOMATIC_MIGRATION);
                geom = get_geometry_from_vehicles(rank, zz, vehicles, SIZE_X, SIZE_Y); //wont work with HSFC
                create_random_sources(rank, -1, SIZE_X, SIZE_Y, geom, &ca_matrix); //wont work with HSFC
                std::fill(window.data_container.begin(), window.data_container.end(), 0);
            }
        }
#endif
        MPI_Barrier(bottom);
        PAR_STOP_TIMING(step_time, bottom);

        if(step % 100 == 0) {
            if(!rank) set_crash_maker_on_out_roads(rcnt > 0, SIZE_X, SIZE_Y, geom, &ca_matrix);
            update_crash_positions(rcnt > 0, SIZE_X, SIZE_Y, crash, &ca_matrix, bottom);
        }

        /****************************************Start printing**************************************************/
        std::vector<Vehicle> all_vehicles;

        if(wsize > 1) {
            gather_elements_on(vehicles, 0, &all_vehicles, datatype.elements_datatype, bottom);
            gather_elements_on(top_vehicles, 0, &all_vehicles, datatype.elements_datatype, bottom);
        } else all_vehicles = vehicles;

        if(!rank) {
            std::cout << "Time for step " << step
                      << "; [TOT " << step_time
                      << ", CPT " << computation_time
                      << ", COM EXCHANGE " << comm_time
                      << ", COM MIGRATE " << migrate_time
                      << "] => " << (100*computation_time/step_time)<<"% CPT "
                      << (100*comm_time/step_time)<<"% COM" << std::endl;

            auto vehicle_matrix_print = to_map(SIZE_X, SIZE_Y, all_vehicles);
            auto img = zzframe( SIZE_X, SIZE_Y, ca_matrix, vehicle_matrix_print);
            auto step_str =  std::to_string(step);
            step_str = std::string(std::to_string(MAX_STEP).length() - step_str.length(), '0') + step_str;
            out.open(step_str+"_waiting_time.txt", std::ofstream::out);
            print_vehicles(out, all_vehicles);
            img.save(( step_str+prefix_fname).c_str());
            out.close();
        }
        // Stop printing
        /********************************************************************************************************/
        step++;
    }

    MPI_Finalize();
    return 0;
}

//
// Created by xetql on 12/11/18.
//

#ifndef NBMPI_UNLOADING_MODEL_HPP
#define NBMPI_UNLOADING_MODEL_HPP

#include <vector>
#include <tuple>
#include <set>
#include <iterator>
#include <sstream>

#include "zoltan_fn.hpp"

namespace load_balancing {
namespace esoteric {

namespace {

template<class A>
A take_out(const size_t data_id, std::vector<A> *data) {
    std::iter_swap(data->begin() + data_id, data->end() - 1);
    A el = *(data->end() - 1);
    data->pop_back();
    return el;
}

template<class InputIterator1, class InputIterator2>
bool contains_at_least_one_of(InputIterator1 first1, InputIterator1 last1,
                              InputIterator2 first2, InputIterator2 last2) {

    while (first1 < last1 && first2 < last2) {
        if (*first1 == *first2)
            return true;
        else if (*first1 < *first2) {
            first1++;
        } else if (*first1 > *first2) {
            first2++;
        }
    }
    return false;
}

}

enum MIGRATION_TYPE { BOTTOM, TOP };
enum MODEL_STATE {
    init,
    started,
    finished,
    ok,
    stopped,
    on_error,
    on_error_too_many_increasing,
    on_error_no_increasing
};

const float SLOPE_THRESHOLD = 2.0f;
const int SLOPE_LISTENING_DURATION = 20;
const float PERCENTAGE_PER_IT = SLOPE_THRESHOLD / SLOPE_LISTENING_DURATION;
const int LB_INITIALIZED = 0, LB_STARTED = 0, LB_FINISHED = 0, LB_OK = 0, LB_ERR_TOO_MANY_INCREASING = 1;
const int COMM_TOP_INDEX = 1, COMM_INCREASING_INDEX = 0;
const int RESET_AFTER = 200;
struct UnloadingModelLocalState {
    mutable MODEL_STATE state = MODEL_STATE::stopped;
    int number_of_elements_within_top = 0,
            number_of_elements_within_bottom = 0,
            number_of_elements_on_loan = 0,
            reset_after = -1;
    bool in_top_partition = true;
    std::vector<int> increasing_cpus;
    int started_at_step = -1;
    float sigma;
    MPI_Comm increasing_comm = MPI_COMM_NULL;

};

/**
 * get communicators for the new LB paradigm
 * @param loading_ranks  all processors must know this.
 * @param world_comm
 * @param top the new communicator for people how will communicate on "top".
 */
// TODO: If everybody has increasing load, then nobody in top, then the execution will fail.
template<class A>
UnloadingModelLocalState
init_unloading_model(int step, float my_load_slope, int my_rank, const std::vector<A> &data, MPI_Comm bottom) {
    // top is set for PE that does not have an increasing load others have undefined.
    // bottom is simply MPI_COMM_WORLD normally
    // easy as fuck right?
    const int TOP_CPU_TAG = 800;
    bool is_not_increasing = my_load_slope < SLOPE_THRESHOLD;
    MODEL_STATE err = MODEL_STATE::init;
    MPI_Comm top, incr = MPI_COMM_NULL;
    MPI_Comm_split(bottom, is_not_increasing ? COMM_TOP_INDEX : COMM_INCREASING_INDEX, my_rank, &top);
    MPI_Group top_gr;
    MPI_Group bottom_gr;
    MPI_Comm_group(bottom, &bottom_gr);
    int bottom_gr_size;
    MPI_Comm_size(bottom, &bottom_gr_size);
    int top_gr_size, incr_gr_size;
    std::vector<int> increasing_cpus;
    if (is_not_increasing) { //top_group
        int top_rank;
        MPI_Comm_rank(top, &top_rank);
        MPI_Comm_group(top, &top_gr);
        MPI_Group_size(top_gr, &top_gr_size);

        if (top_gr_size == bottom_gr_size) {
            err = MODEL_STATE::on_error_no_increasing;
        } else {
            std::vector<int> top_ranks(top_gr_size);
            std::iota(top_ranks.begin(), top_ranks.end(), 0);
            MPI_Group_translate_ranks(top_gr, top_ranks.size(), &top_ranks.front(), bottom_gr, &top_ranks.front());

            MPI_Group increasing_gr;
            MPI_Group_difference(bottom_gr, top_gr, &increasing_gr);
            int increasing_gr_size;
            MPI_Group_size(increasing_gr, &increasing_gr_size);

            if (increasing_gr_size >= bottom_gr_size / 2)
                err = MODEL_STATE::on_error_too_many_increasing;
            else {
                increasing_cpus.resize(increasing_gr_size);
                std::iota(increasing_cpus.begin(), increasing_cpus.end(), 0);
                MPI_Group_translate_ranks(increasing_gr, increasing_cpus.size(), &increasing_cpus.front(), bottom_gr,
                                          &increasing_cpus.front());
            }
        }
    } else { //increasing_group
        MPI_Comm_size(top, &incr_gr_size);
        if (incr_gr_size >= bottom_gr_size / 2) {
            err = MODEL_STATE::on_error_too_many_increasing;
        } else {
            MPI_Comm_group(top, &top_gr);
            increasing_cpus.resize(incr_gr_size);
            std::iota(increasing_cpus.begin(), increasing_cpus.end(), 0);
            MPI_Group_translate_ranks(top_gr, increasing_cpus.size(), &increasing_cpus.front(), bottom_gr,
                                      &increasing_cpus.front());
            top_gr_size = bottom_gr_size - incr_gr_size;
            MPI_Group incr_gr;
            MPI_Group_incl(bottom_gr, increasing_cpus.size(), &increasing_cpus.front(), &incr_gr);
            MPI_Comm_create_group(bottom, incr_gr, 0, &incr);
        }
    }

    std::sort(increasing_cpus.begin(), increasing_cpus.end());
    float sigma = -1;
    if (err != MODEL_STATE::init) {
        increasing_cpus.clear();
    } else {
        int total_nb_vehicles, my_nb_vehicles = data.size();
        MPI_Allreduce(&my_nb_vehicles, &total_nb_vehicles, 1, MPI_INT, MPI_SUM, bottom);
        int mu = total_nb_vehicles / bottom_gr_size;
        float max_slope;
        MPI_Allreduce(&my_load_slope, &max_slope, 1, MPI_FLOAT, MPI_MAX, bottom);
        sigma = mu / max_slope;
    }
    return {err, top_gr_size, bottom_gr_size, 0, -1, is_not_increasing, increasing_cpus, step, sigma,
            is_not_increasing ? MPI_COMM_NULL : incr};
}


template<class A>
Zoltan_Struct* start_unloading_model(std::vector<A> *data_bottom, // becomes bottom
                                     std::vector<A> *data_top,
                                     const UnloadingModelLocalState &model_state,
                                     const MPI_Datatype datatype,
                                     MPI_Comm bottom) {
    const int TAG = 900;
    std::vector<A> top_mesh_data;
    const std::vector<int> &increasing_cpus = model_state.increasing_cpus;
    Zoltan_Struct *zz_top = nullptr;
    if(increasing_cpus.size() == 0) return zz_top;
    switch (model_state.state) {
        case MODEL_STATE::init:
            break;
        default:
            return zz_top;
    }

    int bottom_rank;
    MPI_Comm_rank(bottom, &bottom_rank);
    int bottom_size;
    MPI_Comm_size(bottom, &bottom_size);

    const bool in_top_partition =
            std::find(increasing_cpus.cbegin(), increasing_cpus.cend(), bottom_rank) == increasing_cpus.cend();

    if (!in_top_partition)
        std::move(data_bottom->begin(), data_bottom->end(), std::back_inserter(top_mesh_data));

    zz_top = zoltan_create_wrapper(ENABLE_AUTOMATIC_MIGRATION, bottom, bottom_size - increasing_cpus.size(), in_top_partition ? 1 : 0);
    tca::zoltan_load_balance(&top_mesh_data, zz_top, ENABLE_AUTOMATIC_MIGRATION);

    if (in_top_partition) {
        std::move(top_mesh_data.begin(), top_mesh_data.end(), std::back_inserter(*data_top));
    } else data_bottom->clear();

    model_state.state = MODEL_STATE::started;

    return zz_top;
}
template<class A>
const std::vector<A> exchange(
        Zoltan_Struct   *zoltan_bottom,
        Zoltan_Struct   *zoltan_top,
        std::vector<A>  *bottom_data,
        std::vector<A>  *top_data, // it is always empty for increasing load cpus
        int* nb_elements_recv, int* nb_elements_sent,
        const std::tuple<int, int, int, int>& geom,
        const std::vector<int>& increasing_cpus,
        const MPI_Datatype datatype,
        MPI_Comm bottom,
        const double cell_size = 1.0) {


    if(zoltan_top == nullptr) {
        return zoltan_exchange_data(zoltan_bottom, bottom_data, nb_elements_recv, nb_elements_sent, datatype, bottom, cell_size);
    }

    *nb_elements_recv = 0;
    *nb_elements_sent = 0;

    int my_bottom_rank; MPI_Comm_rank(bottom, &my_bottom_rank);
    int wsize; MPI_Comm_size(bottom, &wsize);

    std::vector<int> num_import_from_procs(wsize);
    std::vector<int> import_from_procs;

    const bool in_top_partition = std::find(increasing_cpus.cbegin(), increasing_cpus.cend(), my_bottom_rank) == increasing_cpus.cend();
    const int EXCHANGE_SEND_TAG_TOP=903, EXCHANGE_SEND_TAG_BOTTOM=904;

    std::vector<A> buffer, remote_data_gathered;

    std::vector<std::vector<A>> data_to_migrate(wsize), data_to_migrate_bottom(wsize);
    std::vector<int> export_gids_b, export_lids_b, export_procs_b, export_gids, export_lids, export_procs;
    std::vector<int> parts; int num_found = 0, num_known = 0;

    num_known = 0;
    size_t data_id = 0;
    std::vector<int> PEs_top(wsize, -1), PEs_bottom(wsize, -1);
    const size_t bot_data_size = bottom_data->size();

    /*********************************************************************************************************
     * COMPUTE DESTINATION FOR DATA ON BOTTOM PROCS
     *********************************************************************************************************/
    std::vector<int> PEs_distinct;
    std::vector<int> PE_attendance(wsize, 0);
    int num_found_proc, num_found_part;
    while (data_id < bot_data_size) {
        Vehicle& v = bottom_data->at(data_id);
        //std::cout << my_bottom_rank << " " << v << std::endl;
        std::pair<double, double> pos_in_double = {(double) v.position.first, (double) v.position.second};

        bool bottom_contain_incr_cpu = !in_top_partition;
        num_found_proc = 0;
        if(inside_the_borders(v, geom, cell_size)){
            /***********************************************************************************
             * Check if a work unit must be shared with someone within the bottom partitioning *
             ***********************************************************************************/
            Zoltan_LB_Box_PP_Assign(zoltan_bottom,
                                    pos_in_double.first - cell_size, pos_in_double.second - cell_size, 0.0,
                                    pos_in_double.first + cell_size, pos_in_double.second + cell_size, 0.0,
                                    &PEs_bottom.front(), &num_found_proc, &parts.front(), &num_found_part);
            bottom_contain_incr_cpu = contains_at_least_one_of(PEs_bottom.begin(), PEs_bottom.begin()+num_found_proc, increasing_cpus.begin(), increasing_cpus.end());
        }
        for (auto iPE = PEs_bottom.begin(); iPE < PEs_bottom.begin()+num_found_proc; iPE++) {
            auto PE = *iPE;
            if(PE >= 0 && !PE_attendance[PE]){
                PE_attendance[PE]++;
                PEs_distinct.push_back(PE);
            }
        }
        if(bottom_contain_incr_cpu) {
            /**********************************************************************************
             * Check if a work unit must be shared with someone within the top partitioning   *
             **********************************************************************************/
            Zoltan_LB_Box_PP_Assign(zoltan_top,
                                    pos_in_double.first - cell_size, pos_in_double.second - cell_size, 0.0,
                                    pos_in_double.first + cell_size, pos_in_double.second + cell_size, 0.0,
                                    &PEs_top.front(), &num_found_proc, &parts.front(), &num_found_part);
            for (auto iPE = PEs_top.begin(); iPE < PEs_top.begin()+num_found_proc; iPE++) {
                auto PE = *iPE;
                if(PE >= 0 && !PE_attendance[PE]){
                    PE_attendance[PE]++;
                    PEs_distinct.push_back(PE);
                }
            }
        }

        /**********************************************************************************
         * Mark data as export for every detected PE                                      *
         **********************************************************************************/
        for(int PE : PEs_distinct) {
            if (PE >= 0 && PE != my_bottom_rank) {
                export_gids.push_back(bottom_data->at(data_id).gid);
                export_lids.push_back(bottom_data->at(data_id).lid);
                export_procs.push_back(PE);
                data_to_migrate.at(PE).push_back(bottom_data->at(data_id));
                num_known++;
            }
        }
        PEs_distinct.clear();
        std::fill(PE_attendance.begin(), PE_attendance.end(), 0);
        data_id++; //if the element must stay with me then check the next one
    }

    /*********************************************************************************************************
     * EXCHANGE DATA WITH BOTTOM
     *********************************************************************************************************/

    // Compute who has to send me something via Zoltan.
    num_import_from_procs.resize(wsize);
    {
        std::vector<MPI_Request> rcv_reqs(wsize);
        std::vector<MPI_Status> statuses(wsize);
        std::vector<int> export_to_procs, size_to_send(wsize, 0);
        std::fill(num_import_from_procs.begin(), num_import_from_procs.end(), 0);
        int i = 0;
        for (int PE = 0; PE < wsize; PE++) {
            MPI_Irecv(&num_import_from_procs[PE], 1, MPI_INT, PE, 1234, bottom, &rcv_reqs[PE]);
            if(!data_to_migrate[PE].empty()) {
                export_to_procs.push_back(PE);
                int send_size = data_to_migrate[PE].size();
                MPI_Send(&send_size, 1, MPI_INT, PE, 1234, bottom);
                i++;
            }
        }
        MPI_Barrier(bottom);
        import_from_procs.clear();
        for (size_t PE = 0; PE < wsize; PE++) {
            int flag; MPI_Status status;
            MPI_Test(&rcv_reqs[PE], &flag, &status);
            if(!flag)
                MPI_Cancel(&rcv_reqs[PE]);
            else {
                import_from_procs.push_back(PE);
            }
        }
    }

    // Compute how many elements I have to import from others, and from whom.

    // if nothing found, nothing to free.

    auto nb_reqs = std::count_if(data_to_migrate.cbegin(), data_to_migrate.cend(), [](auto buf){return !buf.empty();});
    int cpt = 0;

    // Send the data to neighbors
    std::vector<MPI_Request> bottom_reqs(nb_reqs);
    for (size_t PE = 0; PE < wsize; PE++) {
        int send_size = data_to_migrate.at(PE).size();
        if (send_size) {
            MPI_Isend(&data_to_migrate.at(PE).front(), send_size, datatype, PE, EXCHANGE_SEND_TAG_BOTTOM, bottom,
                      &bottom_reqs[cpt]);
            cpt++;
        }
    }

    // Import the data from neighbors
    std::unordered_map<int, std::vector<A>> bottom_data_tmp;
    for (int proc_id : import_from_procs) {
        auto size = num_import_from_procs[proc_id];
        bottom_data_tmp[proc_id].resize(size);
        MPI_Recv(&bottom_data_tmp[proc_id].front(), size, datatype, proc_id, EXCHANGE_SEND_TAG_BOTTOM, bottom, MPI_STATUS_IGNORE);
    }

    /*********************************************************************************************************
     * COMPUTE DESTINATION FOR TOP DATA ON TOP PROCS
     *********************************************************************************************************/

    std::fill(PEs_top.begin(), PEs_top.end(), -1);
    std::fill(PEs_bottom.begin(), PEs_bottom.end(), -1);
    std::fill(data_to_migrate.begin(), data_to_migrate.end(), std::vector<Vehicle>() );
    std::fill(num_import_from_procs.begin(), num_import_from_procs.end(), 0);
    export_gids.clear();
    export_procs.clear();
    export_lids.clear();
    num_known = 0;
    data_id = 0;
    const size_t top_data_size = top_data->size();

    if(in_top_partition) { //can't be executed by bottom-only PEs
        while (data_id < top_data_size) {
            int num_found_proc, num_found_part;
            Vehicle& v = top_data->at(data_id); //!\\ TOP DATA HERE
            /**********************************************************************************
             * Check if a work unit must be shared with someone within the bottom partitioning*
             **********************************************************************************/
            std::pair<double, double> pos_in_double = {(double) v.position.first, (double) v.position.second};
            Zoltan_LB_Box_PP_Assign(zoltan_bottom,
                                    pos_in_double.first - cell_size, pos_in_double.second - cell_size, 0.0,
                                    pos_in_double.first + cell_size, pos_in_double.second + cell_size, 0.0,
                                    &PEs_bottom.front(), &num_found_proc, &parts.front(), &num_found_part);

            // Erase all the PEs that did not send a data that interact with the current vehicle
            std::vector<int> PEs_bottom2(num_found_proc, -1);
            for(int i = 0; i < num_found_proc; ++i){
                int PE = PEs_bottom[i];
                if(exists(bottom_data_tmp, PE)) {
                    bool must_send = std::none_of(bottom_data_tmp[PE].cbegin(), bottom_data_tmp[PE].cend(),[&v, &cell_size](const auto &e) {
                        return distance2(e, v) <= cell_size;
                    });
                    if(must_send) PEs_bottom2[i] = PE;
                }
            }

            PEs_distinct.clear();
            std::fill(PE_attendance.begin(), PE_attendance.end(), 0);
            for (auto iPE = PEs_bottom2.begin(); iPE < PEs_bottom2.end(); iPE++) {
                auto PE = *iPE;
                if(PE >= 0 && !PE_attendance[PE]) {
                    PE_attendance[PE]++;
                    PEs_distinct.push_back(PE);
                }
            }

            /**********************************************************************************
             * Check if a work unit must be shared with someone within the top partitioning   *
             **********************************************************************************/
            Zoltan_LB_Box_PP_Assign(zoltan_top,
                                    pos_in_double.first - cell_size, pos_in_double.second - cell_size, 0.0,
                                    pos_in_double.first + cell_size, pos_in_double.second + cell_size, 0.0,
                                    &PEs_top.front(), &num_found_proc, &parts.front(), &num_found_part);
            for (auto iPE = PEs_top.begin(); iPE < PEs_top.begin()+num_found_proc; iPE++) {
                auto PE = *iPE;
                if(PE >= 0 && !PE_attendance[PE]){
                    PE_attendance[PE]++;
                    PEs_distinct.push_back(PE);
                }
            }

            /**********************************************************************************
             * Mark data as export for every detected PE                                      *
             **********************************************************************************/
            for (int PE : PEs_distinct) {
                if (PE >= 0 && PE != my_bottom_rank) {
                    export_gids.push_back(top_data->at(data_id).gid);
                    export_lids.push_back(top_data->at(data_id).lid);
                    export_procs.push_back(PE);
                    data_to_migrate.at(PE).push_back(top_data->at(data_id));
                    num_known++;
                }
            }
            data_id++; //if the element must stay with me then check the next one
        }
    }

    /*********************************************************************************************************
     * EXCHANGE DATA TOP BOTTOM
     *********************************************************************************************************/

    import_from_procs.clear();
    num_import_from_procs.resize(wsize);
    {
        std::vector<MPI_Request> rcv_reqs(wsize);
        std::vector<MPI_Status> statuses(wsize);
        std::fill(num_import_from_procs.begin(), num_import_from_procs.end(), 0);
        int i = 0;
        for (int PE = 0; PE < wsize; PE++) {
            MPI_Irecv(&num_import_from_procs[PE], 1, MPI_INT, PE, 4321, bottom, &rcv_reqs[PE]);
            if(!data_to_migrate[PE].empty()) {
                int send_size = data_to_migrate[PE].size();
                MPI_Send(&send_size, 1, MPI_INT, PE, 4321, bottom);
                i++;
            }
        }
        MPI_Barrier(bottom);
        import_from_procs.clear();
        for (size_t PE = 0; PE < wsize; PE++) {
            int flag; MPI_Status status;
            MPI_Test(&rcv_reqs[PE], &flag, &status);
            if(!flag)
                MPI_Cancel(&rcv_reqs[PE]);
            else {
                import_from_procs.push_back(PE);
            }
        }
    }

    nb_reqs = std::count_if(data_to_migrate.cbegin(), data_to_migrate.cend(), [](auto buf){return !buf.empty();});
    cpt = 0;

    // Send the data to neighbors
    std::vector<MPI_Request> top_reqs(nb_reqs);
    for (size_t PE = 0; PE < wsize; PE++) {
        int send_size = data_to_migrate.at(PE).size();
        if (send_size) {
            MPI_Isend(&data_to_migrate.at(PE).front(), send_size, datatype, PE, EXCHANGE_SEND_TAG_TOP, bottom,
                      &top_reqs[cpt]);
            cpt++;
        }
    }

    // Import the data from neighbors
    for (int proc_id : import_from_procs) {
        auto size = num_import_from_procs[proc_id];
        buffer.resize(size);
        MPI_Recv(&buffer.front(), size, datatype, proc_id, EXCHANGE_SEND_TAG_TOP, bottom, MPI_STATUS_IGNORE);
        std::move(buffer.begin(), buffer.end(), std::back_inserter(remote_data_gathered));
    }

    for(auto& data : bottom_data_tmp) {
        std::move(data.second.begin(), data.second.end(), std::back_inserter(remote_data_gathered));
    }

    /// Wait on my requests to complete
    MPI_Waitall(bottom_reqs.size(), &bottom_reqs.front(), MPI_STATUSES_IGNORE);
    MPI_Waitall(top_reqs.size(),    &top_reqs.front(),    MPI_STATUSES_IGNORE);

    return remote_data_gathered;
}

template<class A>
void migrate(Zoltan_Struct *zoltan_bottom,
             Zoltan_Struct *zoltan_top, // both in the same comm
             std::vector<A> *bottom_data,
             std::vector<A> *top_data, // it is always empty for increasing load cpus
             const std::vector<int> &increasing_cpus,
             const MPI_Datatype datatype,
             MPI_Comm bottom) {

    if (zoltan_top == nullptr) {
        zoltan_migrate_particles(zoltan_bottom, bottom_data, datatype, bottom);
        return;
    }

    int wsize;
    MPI_Comm_size(bottom, &wsize);
    int my_bottom_rank;
    MPI_Comm_rank(bottom, &my_bottom_rank);

    const int BOTTOM_SEND_TAG = 901, TOP_SEND_TAG = 902;
    const bool in_top_partition =
            std::find(increasing_cpus.cbegin(), increasing_cpus.cend(), my_bottom_rank) == increasing_cpus.cend();

    std::vector<std::vector<A>> data_to_migrate_top(wsize), data_to_migrate_bottom(wsize);
    std::vector<int> export_gids_b, export_lids_b, export_procs_b,
            export_gids_t, export_lids_t, export_procs_t;

    size_t data_id = 0;
    int PE, part, num_known_top = 0, num_known_bottom = 0;

    // Computing destination for bottom particles
    while (data_id < bottom_data->size()) {
        std::vector<double> pos_in_double = {(double) bottom_data->at(data_id).position.first,
                                             (double) bottom_data->at(data_id).position.second};

        Zoltan_LB_Point_PP_Assign(zoltan_bottom, &pos_in_double.front(), &PE, &part);
        if (PE != my_bottom_rank) {
            export_gids_b.push_back(bottom_data->at(data_id).gid);
            export_lids_b.push_back(bottom_data->at(data_id).lid);
            export_procs_b.push_back(PE);
            data_to_migrate_bottom.at(PE).push_back(take_out(data_id, bottom_data));
            num_known_bottom++;
        } else data_id++; //if the element must stay with me then check the next one
    }

    int my_top_rank = -1;
    // I have top data to migrate because I am not an "increasing" cpu
    if (in_top_partition) {
        data_id = 0;
        // Computing destination for top particles
        while (data_id < top_data->size()) {
            std::vector<double> pos_in_double = {(double) top_data->at(data_id).position.first,
                                                 (double) top_data->at(data_id).position.second};
            Zoltan_LB_Point_PP_Assign(zoltan_bottom, &pos_in_double.front(), &PE, &part);
            if (PE != my_bottom_rank) {
                if (std::find(increasing_cpus.cbegin(), increasing_cpus.cend(), PE) != increasing_cpus.cend()) {
                    Zoltan_LB_Point_PP_Assign(zoltan_top, &pos_in_double.front(), &PE, &part);
                    if (PE != my_bottom_rank) {
                        export_gids_t.push_back(top_data->at(data_id).gid);
                        export_lids_t.push_back(top_data->at(data_id).lid);
                        export_procs_t.push_back(PE);
                        data_to_migrate_top.at(PE).push_back(take_out(data_id, top_data));
                        num_known_top++;
                    } else data_id++;
                } else {
                    export_gids_b.push_back(top_data->at(data_id).gid);
                    export_lids_b.push_back(top_data->at(data_id).lid);
                    export_procs_b.push_back(PE);
                    data_to_migrate_bottom.at(PE).push_back(take_out(data_id, top_data));
                    num_known_bottom++;
                }
            } else data_id++; // if the element must stay with me then check the next one
        }
    }

    /*************************************************/
    /** Variables for send and recv                 **/
    /*************************************************/

    ZOLTAN_ID_PTR found_gids, found_lids;
    int *found_procs, *found_parts, num_found = 0, ierr;
    std::vector<int> num_import_from_procs_t(wsize), num_import_from_procs_b(wsize);
    std::vector<int> import_from_procs_t, import_from_procs_b;

    /*************************************************/
    /** Sending part to BOTTOM                      **/
    /*************************************************/

    auto known_gids = (ZOLTAN_ID_PTR) &export_gids_b.front();
    auto known_lids = (ZOLTAN_ID_PTR) &export_lids_b.front();

    ierr = Zoltan_Invert_Lists(zoltan_bottom, num_known_bottom, known_gids, known_lids, &export_procs_b[0],
                               &export_procs_b[0],
                               &num_found, &found_gids, &found_lids, &found_procs, &found_parts);

    for (size_t i = 0; i < num_found; ++i) {
        num_import_from_procs_b[found_procs[i]]++;
        if (std::find(import_from_procs_b.begin(), import_from_procs_b.end(), found_procs[i]) ==
            import_from_procs_b.end())
            import_from_procs_b.push_back(found_procs[i]);
    }

    if (num_found > 0)
        Zoltan_LB_Free_Part(&found_gids, &found_lids, &found_procs, &found_parts);

    int nb_reqs = 0;
    for (auto buf: data_to_migrate_bottom)
        if (!buf.empty()) nb_reqs++;

    int cpt = 0;
    std::vector<MPI_Request> reqs_b(nb_reqs);
    for (size_t PE = 0; PE < wsize; PE++) {
        int send_size = data_to_migrate_bottom.at(PE).size();
        if (send_size) {
            MPI_Isend(&data_to_migrate_bottom.at(PE).front(), send_size, datatype, PE,
                      BOTTOM_SEND_TAG, bottom, &reqs_b[cpt]);

            cpt++;
        }
    }

    /*************************************************/
    /** Sending part, to TOP                         */
    /*************************************************/

    //Actually, I take part of the rest of the computation, so I may have data to send to the top comm
    std::vector<MPI_Request> reqs_t;

    known_gids = (ZOLTAN_ID_PTR) &export_gids_t.front();
    known_lids = (ZOLTAN_ID_PTR) &export_lids_t.front();

    ierr = Zoltan_Invert_Lists(zoltan_top, num_known_top, known_gids, known_lids, &export_procs_t[0],
                               &export_procs_t[0],
                               &num_found, &found_gids, &found_lids, &found_procs, &found_parts);

    if (in_top_partition) {

        for (size_t i = 0; i < num_found; ++i) {
            num_import_from_procs_t[found_procs[i]]++;
            if (std::find(import_from_procs_t.begin(), import_from_procs_t.end(), found_procs[i]) ==
                import_from_procs_t.end())
                import_from_procs_t.push_back(found_procs[i]);
        }

        if (num_found > 0)
            Zoltan_LB_Free_Part(&found_gids, &found_lids, &found_procs, &found_parts);

        nb_reqs = 0;
        for (auto buf: data_to_migrate_top) {
            if (!buf.empty()) nb_reqs++;
        }

        cpt = 0;
        reqs_t.resize(nb_reqs);
        for (size_t PE = 0; PE < wsize; PE++) {
            int send_size = data_to_migrate_top.at(PE).size();
            if (send_size) {
                MPI_Isend(&data_to_migrate_top.at(PE).front(), send_size, datatype, PE,
                          TOP_SEND_TAG, bottom, &reqs_t[cpt]);
                cpt++;
            }
        }
    }

    /*************************************************/
    /** Receiving part from BOTTOM                  **/
    /*************************************************/

    std::vector<A> buffer;
    for (int proc_id : import_from_procs_b) {
        size_t size = num_import_from_procs_b[proc_id];
        buffer.resize(size);
        MPI_Recv(&buffer.front(), size, datatype, proc_id, BOTTOM_SEND_TAG, bottom, MPI_STATUS_IGNORE);
        std::move(buffer.begin(), buffer.end(), std::back_inserter(*bottom_data));
    }

    /*************************************************/
    /** Receiving part from TOP                      */
    /*************************************************/
    //I don't have an increasing load so I may receive some data from other "top" cpu.
    if (in_top_partition) {
        for (int proc_id : import_from_procs_t) {
            size_t size = num_import_from_procs_t[proc_id];
            buffer.resize(size);
            MPI_Recv(&buffer.front(), size, datatype, proc_id, TOP_SEND_TAG, bottom, MPI_STATUS_IGNORE);
            std::move(buffer.begin(), buffer.end(), std::back_inserter(*top_data));
        }
    }

    // Waiting on requests
    MPI_Waitall(reqs_b.size(), &reqs_b.front(), MPI_STATUSES_IGNORE);
    if (in_top_partition)
        MPI_Waitall(reqs_t.size(), &reqs_t.front(), MPI_STATUSES_IGNORE);

    // Update local ids
    size_t i = 0;
    const int nb_data_b = bottom_data->size();
    for (auto &e : *bottom_data) {
        e.lid = i;
        i++;
    }
    i = 0;
    for (auto &e: *top_data) {
        e.lid = i + nb_data_b;
        i++;
    }

}

/**
 * Stop the unloading model, returning back to normal state with a new partitioning scheme after sigma iterations.
 * Zoltan_top is destroyed and memory freed.
 * @tparam A
 * @param current_step
 * @param zoltan_bottom
 * @param zoltan_top
 * @param bottom_data
 * @param top_data
 * @param model_state
 * @param bottom
 */
template<class A>
bool stop_unloading_model(int current_step,
                          Zoltan_Struct *zoltan_bottom,
                          Zoltan_Struct *zoltan_top, // both are in the same comm
                          std::vector<A> *bottom_data,
                          std::vector<A> *top_data, // it is always empty for increasing load cpus
                          const UnloadingModelLocalState &model_state, MPI_Comm bottom) {
    /** check if the state of the system meets the requirements for resetting **/
    /*if(!model_state.in_top_partition) { //local deletion
        float ratio = bottom_data->size()/model_state.number_of_elements_on_loan;
        std::vector<float> allratios(model_state.increasing_cpus.size());
        MPI_Allgather(&ratio, 1, MPI_FLOAT, &allratios, 1, MPI_FLOAT, model_state.increasing_comm);
        int size, rank;
        MPI_Comm_size(bottom,&size); MPI_Comm_rank(bottom, &rank);
        if (std::all_of(allratios.cbegin(), allratios.cend(), [](auto v){return v > 0.7;})) { // if I recovered 60% of my original state
            auto lower = std::lower_bound(model_state.increasing_cpus.begin(), model_state.increasing_cpus.end(), rank);
            auto idx   = std::distance(model_state.increasing_cpus.begin(), lower);
        }
    } else {
    }
    */
    int r;
    switch (model_state.state) {
        case MODEL_STATE::started:
            if ((model_state.started_at_step + model_state.sigma) <= current_step) { //global deletion
                MPI_Comm_rank(bottom, &r);
                Zoltan_Destroy(&zoltan_top);
                std::copy(top_data->begin(), top_data->end(), std::back_inserter(*bottom_data));
                top_data->clear();
                tca::zoltan_load_balance(bottom_data, zoltan_bottom, ENABLE_AUTOMATIC_MIGRATION);
                model_state.state = MODEL_STATE::finished;
                return true;
            }
        default:
            return false;
    }
}

} // end of namespace esoteric
} // end of namespace load_balancing

#endif //NBMPI_UNLOADING_MODEL_HPP

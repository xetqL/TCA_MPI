//
// Created by xetql on 18.12.17.
//

#ifndef NBMPI_LOADBALANCER_HPP
#define NBMPI_LOADBALANCER_HPP

#include <map>
#include <mpi.h>

#include "spatial_elements.hpp"
#include "zoltan_fn.hpp"

namespace load_balancing {

    template<int N>
    inline void gather_elements_on(const int world_size,
                                   const int my_rank,
                                   const int nb_elements,
                                   const std::vector<elements::Element<N>> &local_el,
                                   const int dest_rank,
                                   std::vector<elements::Element<N>> &dest_el,
                                   const MPI_Datatype &sendtype,
                                   const MPI_Comm &comm) {
        int nlocal = local_el.size();
        std::vector<int> counts(world_size, 0), displs(world_size, 0);
        MPI_Gather(&nlocal, 1, MPI_INT, &counts.front(), 1, MPI_INT, dest_rank, comm);
        for (int cpt = 0; cpt < world_size; ++cpt) displs[cpt] = cpt == 0 ? 0 : displs[cpt - 1] + counts[cpt - 1];
        if (my_rank == dest_rank) dest_el.resize(nb_elements);
        MPI_Gatherv(&local_el.front(), nlocal, sendtype,
                    &dest_el.front(), &counts.front(), &displs.front(), sendtype, dest_rank, comm);
    }

    using ProcessingElementID=int;

    namespace geometric {

        template<int N>
        const std::vector<elements::Element<N>> __zoltan_exchange_data(const std::vector<elements::Element<N>> &data,
                                                                     Zoltan_Struct *load_balancer,
                                                                     const partitioning::CommunicationDatatype datatype,
                                                                     const MPI_Comm LB_COMM,
                                                                     int &nb_elements_recv,
                                                                     int &nb_elements_sent,
                                                                     double cell_size = 0.007) {
            int wsize;
            MPI_Comm_size(LB_COMM, &wsize);
            int caller_rank;
            MPI_Comm_rank(LB_COMM, &caller_rank);

            std::vector<elements::Element<N>> buffer;
            std::vector<elements::Element<N>> remote_data_gathered;
            // Get the neighbors
            std::vector<std::vector<elements::Element<N>>> data_to_migrate(wsize);

            //check within the remaining elements which belong to the current PE
            size_t data_id = 0;
            while (data_id < data.size()) {
                auto el = data.at(data_id).position;
                std::vector<int> procs;
                int nprocs;
                if (N == 3)
                    Zoltan_LB_Box_Assign(load_balancer,
                                         el.at(0) - cell_size, el.at(1) - cell_size, el.at(2) - cell_size,
                                         el.at(0) + cell_size, el.at(1) + cell_size, el.at(2) + cell_size,
                                         &procs.front(), &nprocs);
                else
                    Zoltan_LB_Box_Assign(load_balancer,
                                         el.at(0) - cell_size, el.at(1) - cell_size, 0,
                                         el.at(0) + cell_size, el.at(1) + cell_size, 0,
                                         &procs.front(), &nprocs);

                for (const int PE : procs) {
                    if (PE == (size_t) caller_rank) continue; //do not check with myself
                    data_to_migrate.at(PE).push_back(
                            data.at(data_id)); // get the value and push it in the "to migrate" vector
                }
                data_id++; //if the element must stay with me then check the next one
            }

            std::vector<MPI_Request> reqs(wsize);
            std::vector<MPI_Status> statuses(wsize);
            int cpt = 0, nb_neighbors = wsize;
            nb_elements_sent = 0;
            for (size_t neighbor_idx = 0; neighbor_idx < wsize; ++neighbor_idx) {   //give all my data to neighbors
                int send_size = data_to_migrate.at(neighbor_idx).size();
                nb_elements_sent += send_size;
                MPI_Isend(&data_to_migrate.at(neighbor_idx).front(), send_size, datatype.elements_datatype,
                          neighbor_idx, 200, LB_COMM, &reqs[cpt]);
                cpt++;
            }
            //Wait that people have sent their particles such that no PE have to know who are its neighbors
            MPI_Barrier(LB_COMM);
            cpt = 0;
            int flag = 1;
            while (flag) {// receive the data in any order
                int source_rank, size;
                MPI_Iprobe(MPI_ANY_SOURCE, 200, LB_COMM, &flag, &statuses[cpt]);
                if (!flag) break;
                source_rank = statuses[cpt].MPI_SOURCE;
                MPI_Get_count(&statuses[cpt], datatype.elements_datatype, &size);
                buffer.resize(size);
                MPI_Recv(&buffer.front(), size, datatype.elements_datatype, source_rank, 200, LB_COMM, &statuses[cpt]);
                std::move(buffer.begin(), buffer.end(), std::back_inserter(remote_data_gathered));
                cpt++;
            }
            //MPI_Waitall(reqs.size(), &reqs.front(), &statuses.front()); //less strict than mpi_barrier
            nb_elements_recv = remote_data_gathered.size();
            return remote_data_gathered;
        }

        template<int N>
        void migrate_zoltan(std::vector<elements::Element<N>> &data, int numImport, int numExport, int *exportProcs,
                            unsigned int *exportGlobalGids,
                            const partitioning::CommunicationDatatype datatype,
                            const MPI_Comm LB_COMM) {
            int wsize;
            MPI_Comm_size(LB_COMM, &wsize);
            int caller_rank;
            MPI_Comm_rank(LB_COMM, &caller_rank);
            std::vector<elements::Element<N> > buffer;
            std::map<int, std::shared_ptr<std::vector<elements::Element<N> > > > data_to_migrate;

            for (int i = 0; i < numExport; ++i)
                if (data_to_migrate.find(exportProcs[i]) == data_to_migrate.end())
                    data_to_migrate[exportProcs[i]] = std::make_shared<std::vector<elements::Element<N>>>();

            for (int i = 0; i < numExport; ++i) {
                auto PE = exportProcs[i];
                auto gid = exportGlobalGids[i];

                //check within the remaining elements which belong to the current PE
                size_t data_id = 0;
                while (data_id < data.size()) {
                    if (gid == (size_t) data[data_id].gid) {
                        //if the current element has to be moved, then swap with the last and pop it out (dont need to move the pointer also)
                        //swap iterator values in constant time
                        std::iter_swap(data.begin() + data_id, data.end() - 1);
                        //get the value and push it in the "to migrate" vector
                        data_to_migrate[PE]->push_back(*(data.end() - 1));
                        //pop the head of the list in constant time
                        data.pop_back();
                    } else data_id++; //if the element must stay with me then check the next one
                }
            }

            std::vector<MPI_Request> reqs(data_to_migrate.size());

            int cpt = 0;
            for (auto const &pe_data : data_to_migrate) {
                int send_size = pe_data.second->size();
                MPI_Isend(&pe_data.second->front(), send_size, datatype.elements_datatype, pe_data.first, 300, LB_COMM,
                          &reqs[cpt]);
                cpt++;
            }
            int collectData = 0;

            while (collectData < numImport) {// receive the data in any order
                int source_rank, size;
                MPI_Status status;
                MPI_Probe(MPI_ANY_SOURCE, 300, LB_COMM, &status);
                source_rank = status.MPI_SOURCE;
                MPI_Get_count(&status, datatype.elements_datatype, &size);
                collectData += size;
                buffer.resize(size);
                MPI_Recv(&buffer.front(), size, datatype.elements_datatype, source_rank, 300, LB_COMM, &status);
                std::move(buffer.begin(), buffer.end(), std::back_inserter(data));

            }
            MPI_Waitall(cpt, &reqs.front(), MPI_STATUSES_IGNORE);

        }

        template<int N>
        const std::vector<elements::Element<N>> zoltan_exchange_data(std::vector<elements::Element<N>> &data,
                                                                     Zoltan_Struct *load_balancer,
                                                                     const partitioning::CommunicationDatatype datatype,
                                                                     const MPI_Comm LB_COMM,
                                                                     int &nb_elements_recv,
                                                                     int &nb_elements_sent,
                                                                     double cell_size = 0.000625) {
            int wsize;
            MPI_Comm_size(LB_COMM, &wsize);
            int caller_rank;
            MPI_Comm_rank(LB_COMM, &caller_rank);

            std::vector<elements::Element<N>> buffer;
            std::vector<elements::Element<N>> remote_data_gathered;

            if(wsize == 1) return remote_data_gathered;

	        std::vector<std::vector<elements::Element<N>>> data_to_migrate(wsize);
            size_t data_id = 0;
            std::vector<int> PEs(wsize, -1);
            int num_found, num_known = 0;
            std::vector<int> export_gids, export_lids, export_procs;

            // so much memory could be allocated here... potentially PE * n * DIM * 44 bytes => so linear in N
            // as DIM << PE <<<< n
            while (data_id < data.size()) {
                auto pos_in_double = functional::map<double>(data.at(data_id).position, [](auto p){return (double) p;});
                Zoltan_LB_Box_Assign(load_balancer,
                                     pos_in_double.at(0) - cell_size,
                                     pos_in_double.at(1) - cell_size,
                                     N == 3 ? pos_in_double.at(2) - cell_size : 0.0,
                                     pos_in_double.at(0) + cell_size,
                                     pos_in_double.at(1) + cell_size,
                                     N == 3 ? pos_in_double.at(2) + cell_size : 0.0,
                                     &PEs.front(), &num_found);

                for(int PE_idx = 0; PE_idx < num_found; PE_idx++) {
                    int PE = PEs[PE_idx];
                    if(PE >= 0) {
                        if (PE != caller_rank) {
                            export_gids.push_back(data.at(data_id).gid);
                            export_lids.push_back(data.at(data_id).lid);
                            export_procs.push_back(PE);
                            //get the value and copy it into the "to migrate" vector
                            data_to_migrate.at(PE).push_back(data.at(data_id));
                            num_known++;
                        }
                    }
                }
                data_id++; //if the element must stay with me then check the next one
            }

            ZOLTAN_ID_PTR known_gids = (ZOLTAN_ID_PTR) &export_gids.front();
            ZOLTAN_ID_PTR known_lids = (ZOLTAN_ID_PTR) &export_lids.front();
            ZOLTAN_ID_PTR found_gids, found_lids;

            int *found_procs, *found_parts;

            // Compute who has to send me something via Zoltan.
            int ierr = Zoltan_Invert_Lists(load_balancer, num_known, known_gids, known_lids, &export_procs[0], &export_procs[0],
                                           &num_found, &found_gids, &found_lids, &found_procs, &found_parts);

            std::vector<int> num_import_from_procs(wsize);
            std::vector<int> import_from_procs;

            // Compute how many elements I have to import from others, and from whom.
            for (size_t i = 0; i < num_found; ++i) {
                num_import_from_procs[found_procs[i]]++;
                if (std::find(import_from_procs.begin(), import_from_procs.end(), found_procs[i]) == import_from_procs.end())
                    import_from_procs.push_back(found_procs[i]);
            }

            // if nothing found, nothing to free.
            if(num_found > 0)
                Zoltan_LB_Free_Part(&found_gids, &found_lids, &found_procs, &found_parts);

            int nb_reqs = 0;
            // Just check how many requests I have to create
            for (auto buf: data_to_migrate) {
                if (!buf.empty()) nb_reqs++;
            }

            int cpt = 0;

            // Send the data to neighbors
            std::vector<MPI_Request> reqs(nb_reqs);
            nb_elements_sent = 0;
            for (size_t PE = 0; PE < wsize; PE++) {
                int send_size = data_to_migrate.at(PE).size();
                if (send_size) {
                    nb_elements_sent += send_size;
                    MPI_Isend(&data_to_migrate.at(PE).front(), send_size, datatype.elements_datatype, PE, 400, LB_COMM,
                              &reqs[cpt]);
                    cpt++;
                }
            }
            // Import the data from neighbors
            nb_elements_recv = 0;
            for (int proc_id : import_from_procs) {
                size_t size = num_import_from_procs[proc_id];
                nb_elements_recv += size;
                buffer.resize(size);
                MPI_Recv(&buffer.front(), size, datatype.elements_datatype, proc_id, 400, LB_COMM, MPI_STATUS_IGNORE);
                std::move(buffer.begin(), buffer.end(), std::back_inserter(remote_data_gathered));
            }

            MPI_Waitall(reqs.size(), &reqs.front(), MPI_STATUSES_IGNORE);
            return remote_data_gathered;
        }
        template<int N>
        void zoltan_migrate_particles(
                std::vector<elements::Element<N>> &data,
                Zoltan_Struct *load_balancer,
                const partitioning::CommunicationDatatype datatype,
                const MPI_Comm LB_COMM) {
            int wsize;
            MPI_Comm_size(LB_COMM, &wsize);
            int caller_rank;
            MPI_Comm_rank(LB_COMM, &caller_rank);

	    if(wsize == 1) return;

            std::vector<std::vector<elements::Element<N>>> data_to_migrate(wsize);

            size_t data_id = 0;
            int PE;
            int num_known = 0;
            std::vector<int> export_gids, export_lids, export_procs;
            while (data_id < data.size()) {
                auto pos_in_double = functional::map<double>(data.at(data_id).position, [](auto p){return (double) p;});
                Zoltan_LB_Point_Assign(load_balancer, &pos_in_double.front(), &PE);
                if (PE != caller_rank) {
                    export_gids.push_back(data.at(data_id).gid);
                    export_lids.push_back(data.at(data_id).lid);
                    export_procs.push_back(PE);
                    //if the current element has to be moved, then swap with the last and pop it out (dont need to move the pointer also)
                    //swap iterator values in constant time
                    std::iter_swap(data.begin() + data_id, data.end() - 1);
                    //get the value and push it in the "to migrate" vector
                    data_to_migrate.at(PE).push_back(*(data.end() - 1));
                    //pop the head of the list in constant time
                    data.pop_back();
                    num_known++;
                } else data_id++; //if the element must stay with me then check the next one
            }

            ZOLTAN_ID_PTR known_gids = (ZOLTAN_ID_PTR) &export_gids.front();
            ZOLTAN_ID_PTR known_lids = (ZOLTAN_ID_PTR) &export_lids.front();
            ZOLTAN_ID_PTR found_gids, found_lids;

            int *found_procs, *found_parts, num_found;

            int ierr = Zoltan_Invert_Lists(load_balancer, num_known, known_gids, known_lids, &export_procs[0], &export_procs[0],
                                &num_found, &found_gids, &found_lids, &found_procs, &found_parts);

            std::vector<int> num_import_from_procs(wsize);
            std::vector<int> import_from_procs;

            for (size_t i = 0; i < num_found; ++i) {
                num_import_from_procs[found_procs[i]]++;
                if (std::find(import_from_procs.begin(), import_from_procs.end(), found_procs[i]) == import_from_procs.end())
                    import_from_procs.push_back(found_procs[i]);
            }

            /* Let's Migrate ma boi ! */

            if(num_found > 0)
                Zoltan_LB_Free_Part(&found_gids, &found_lids, &found_procs, &found_parts);

            int nb_reqs = 0;
            for (auto buf: data_to_migrate) {
                if (!buf.empty()) nb_reqs++;
            }

            int cpt = 0;
            std::vector<MPI_Request> reqs(nb_reqs);
            for (size_t PE = 0; PE < wsize; PE++) {
                int send_size = data_to_migrate.at(PE).size();
                if (send_size) {
                    MPI_Isend(&data_to_migrate.at(PE).front(), send_size, datatype.elements_datatype, PE, 300, LB_COMM,
                              &reqs[cpt]);
                    cpt++;
                }
            }
            std::vector<elements::Element<N>> buffer;
            for (int proc_id : import_from_procs) {
                size_t size = num_import_from_procs[proc_id];
                buffer.resize(size);
                MPI_Recv(&buffer.front(), size, datatype.elements_datatype, proc_id, 300, LB_COMM, MPI_STATUS_IGNORE);
                std::move(buffer.begin(), buffer.end(), std::back_inserter(data));
            }

            const int nb_data = data.size();
            for(int i = 0; i < nb_data; ++i) data[i].lid = i;

            MPI_Waitall(reqs.size(), &reqs.front(), MPI_STATUSES_IGNORE);

        }
    }
}

#endif //NBMPI_LOADBALANCER_HPP
